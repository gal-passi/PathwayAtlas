import pandas as pd
from utils import KeggApi, multiprocess_task
from definitions import *
from os.path import join as pjoin
from utils import save_obj, load_obj, read_in_chunks
import os
from functools import partial


def create_kegg_gene(gene_id_set):
    gene_id = next(iter(gene_id_set))  # since each task is a set with one gene_id
    return KeggGene(gene_id)


class KeggNetwork:
    """Object to represent Kegg Modules or Pathways"""

    def __init__(self, kegg_id, network_type):
        """Constructor for KeggPathway"""
        self.id, self.type, self.gene_list, self._len = None, None, None, {'aa': 0, 'na': 0}
        self._directory = pjoin(KEGG_PATHWAY_OBJECTS_PATH, kegg_id + '.pickle')
        if os.path.exists(self._directory):
            load_obj(self, self._directory, name=kegg_id)
        else:
            self._create_new_instance(kegg_id, network_type)

    def __len__(self):
        return self._len['na']

    def _create_new_instance(self, kegg_id, network_type):
        assert network_type.lower() in NETWORK_TYPES, NETWORK_TYPE_ERROR
        kegg_api = KeggApi()
        self.id, self.type = kegg_id, network_type.lower()
        self.gene_list = kegg_api.get_gene_list(kegg_id)
        gene_objects = []

        def _callback(collector, kegg_gene):
            collector.append(kegg_gene)

        callback = partial(_callback, gene_objects)
        tasks = [(gene_id, ) for gene_id in self.gene_list]  # tasks must be iterable of iterables
        #  too many workers may overload Kegg servers
        multiprocess_task(tasks=tasks, target=create_kegg_gene, callback=callback, workers=6)
        for gene in gene_objects:
            self._len['aa'] += gene.length('aa')
            self._len['na'] += gene.length('na')
        save_obj(self, self._directory)

    @property
    def genes(self):
        """
        iterator for pathway gene objects
        :return: iteratoe KeggGene
        """
        for gene in self.gene_list:
            yield KeggGene(gene)

    def all_snvs(self, outpath='', index=True):
        """
        creates a DataFrame of all single nucleotide variants in the path
        :param index: bool include index in DataFrame
        :param outpath: if not given will be saved to default path at KEGG_PATHWAY_MUTATIONS_PATH
        :return: pandas DataFrame
        """

        tasks = [{gene} for gene in self.genes]
        target = lambda gene: gene.all_snvs()
        collector = []

        def callback_(array, data):
            #print(f"called! {len(array)}")
            array.append(data)

        callback = partial(callback_, collector)
        multiprocess_task(tasks=tasks, target=target, callback=callback)
        all_snvs = pd.concat(collector)
        if not outpath:
            outpath = pjoin(KEGG_PATHWAY_MUTATIONS_PATH, f"{self.id}.csv")
        all_snvs.to_csv(outpath, index=index)
        return all_snvs


class KeggGene:
    """Gene instance for KEGG pathway"""

    def __init__(self, kegg_id, default_init=True):
        """Constructor for Protein"""
        self.kegg_id, self.uniprot_id, self.ref_names = kegg_id, None, None
        self.na_seq, self.aa_seq, self.chr, self.start, self.end = None, None, None, None, None
        self.coding_type = None
        self._dir_name = kegg_id.replace(':', '_')
        self._directory = pjoin(KEGG_GENES_PATH, self._dir_name + '.pickle')
        if os.path.exists(self._directory):
            load_obj(self, self._directory, name=kegg_id)
        elif default_init:
            self._create_new_instance(kegg_id)

    def __len__(self):
        return len(self.aa_seq)

    def _create_new_instance(self, kegg_id):
        kegg_api = KeggApi()
        self.uniprot_id = kegg_api.convert_gene_names(kegg_id)  # dict
        self.aa_seq = kegg_api.gene_seq(kegg_id, 'aaseq')[kegg_id]
        self.na_seq = kegg_api.gene_seq(kegg_id, 'ntseq')[kegg_id]
        save_obj(self, self._directory)

    def create_from_dict(self, data):
        self.__dict__ = data
        self._dir_name = self.kegg_id.replace(':', '_')
        self._directory = pjoin(KEGG_GENES_PATH, self._dir_name + '.pickle')
        save_obj(self, self._directory)

    @property
    def uid(self):
        """
        primary uniprot id
        :return: str
        """
        if self.uniprot_id and 'primary' in self.uniprot_id:
            return self.uniprot_id['primary']
        return None

    @property
    def alias_uid(self):
        """
        alias uniprot ids does not return the main id
        :return: set
        """
        if self.uniprot_id and 'secondary' in self.uniprot_id:
            return set(self.uniprot_id['secondary'])
        return set()

    def length(self, seq='aa'):
        """
        :param seq: aa | na
        :return: return length of amino acid or nucleic acid sequence
        """
        return len(self.aa_seq) if seq == 'aa' else len(self.na_seq)

    def all_snvs(self, outpath='', index=False):
        """
        Creates a DataFrame of all nonsynonymous single nucleotide variants (SNVs) in the gene.

        Notes:
            - The last codon (3 nucleotides) is skipped (assumed to be the stop codon).
            - Any leftover 1 or 2 nucleotides at the end are ignored — only full codons are processed.

        :param index: bool, include index column in DataFrame if True
        :param outpath: str, if provided, saves the DataFrame as a CSV to this path
        :return: pandas DataFrame with SNV data
        """

        def read_in_chunks(seq, chunk_size=3):
            """Yield only full codons (3 bases)."""
            for i in range(0, len(seq) - (len(seq) % chunk_size), chunk_size):
                yield seq[i:i + chunk_size]

        row_data = lambda idx, ref_na, alt_na, ref_aa, alt_aa: \
            ['-', idx, idx, ref_na, alt_na, self.uid, f'{ref_aa}{idx}{alt_aa}']

        mutate_codon = lambda codon, idx, alt: codon[:idx] + alt + codon[idx + 1:]

        df = pd.DataFrame(columns=FAMANALYSIS_COLUMNS)

        for chunk, codon in enumerate(read_in_chunks(self.na_seq[:-3], chunk_size=CODON_LENGTH)):  # skip stop codon
            codon = codon.lower()  # ensure lowercase for consistency with CODON_TRANSLATOR
            if codon not in CODON_TRANSLATOR:
                continue  # skip unknown or invalid codons
            ref_aa = CODON_TRANSLATOR[codon]

            for idx, ref_na in enumerate(codon):
                if ref_na not in NA_CHANGE:
                    continue  # skip invalid nucleotide
                index = (CODON_LENGTH * chunk) + idx

                for alt_na in NA_CHANGE[ref_na]:
                    alt_codon = mutate_codon(codon, idx, alt_na)
                    alt_codon = alt_codon.lower()
                    if alt_codon not in CODON_TRANSLATOR:
                        continue  # skip invalid mutated codons
                    alt_aa = CODON_TRANSLATOR[alt_codon]
                    if alt_aa == ref_aa or alt_aa == STOP_AA:                       #  IS THE SECOND PART NEEDED?
                        continue  # ignore synonymous and nonsense variants

                    df.loc[len(df)] = row_data(index, ref_na, alt_na, ref_aa, alt_aa)

        if outpath:
            df.to_csv(outpath, index=index)

        return df
