import pandas as pd
from utils import KeggApi
from definitions import *
from os.path import join as pjoin
from utils import save_obj, load_obj
import os
from tqdm import tqdm

def gene_snvs_wrapper(gene):
    return gene.all_snvs()

def create_kegg_gene(gene_id_set):
    gene_id = next(iter(gene_id_set))  # since each task is a set with one gene_id
    return KeggGene(gene_id)


class KeggNetwork:
    """Object to represent KEGG Modules or Pathways using precomputed gene SNV files."""

    def __init__(self, kegg_id, network_type):
        self.id, self.type = kegg_id, network_type.lower()
        self._dict_path = pjoin(KEGG_PATHWAY_OBJECTS_PATH, f"{self.id}.pickle")

        assert self.type in NETWORK_TYPES, NETWORK_TYPE_ERROR

        if not os.path.exists(self._dict_path):
            raise FileNotFoundError(f"Missing gene SNV mapping for pathway/module: {self.id}")

        self.gene_snv_map: dict = pd.read_pickle(self._dict_path)  # {kegg_id: path_to_snv_csv}
        self.gene_list = list(self.gene_snv_map.keys())

    def __len__(self):
        return len(self.gene_list)

    @property
    def genes(self):
        """Yield KeggGene instances for all genes in the network."""
        for gene_id in self.gene_list:
            yield KeggGene(gene_id)

    def all_snvs(self, outpath='', index=True):
        """
        Load precomputed SNVs for all genes in the pathway.
        :param index: bool, include index in final CSV
        :param outpath: optional output CSV path (defaults to KEGG_PATHWAY_MUTATIONS_PATH/<id>.csv)
        :return: DataFrame of concatenated SNVs
        """
        collector = []
        for gene_id in tqdm(self.gene_list, desc=f"Reading SNVs for {self.id}", unit="gene"):
            snv_file = self.gene_snv_map.get(gene_id)
            if snv_file and os.path.exists(snv_file):
                try:
                    df = pd.read_csv(snv_file)
                    if not df.empty:
                        collector.append(df)
                except Exception as e:
                    print(f"[Error] Reading {snv_file}: {e}")
            else:
                print(f"[Warning] SNV file not found for {gene_id}")

        if not collector:
            print(f"[Warning] No SNVs found for pathway {self.id}")
            return pd.DataFrame(columns=FAMANALYSIS_COLUMNS)

        all_snvs = pd.concat(collector, ignore_index=True)

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
            - If len(self.na_seq) % 3 != 0, we return empty DataFrame.

        :param index: bool, include index column in DataFrame if True
        :param outpath: str, if provided, saves the DataFrame as a CSV to this path
        :return: pandas DataFrame with SNV data
        """
        # SKIP if the nucleic acid sequence is not a multiple of 3
        if len(self.na_seq) % 3 != 0:
            print(f"Gene: {self.kegg_id} has {len(self.na_seq)} nucleotides, <!%3==0>!")
            return pd.DataFrame(columns=FAMANALYSIS_COLUMNS)

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
                    if alt_aa == ref_aa:            # we include stop codon mutation in the meanwhile
                        continue  # ignore synonymous and nonsense variants

                    df.loc[len(df)] = row_data(index, ref_na, alt_na, ref_aa, alt_aa)

        if outpath:
            df.to_csv(outpath, index=index)

        return df
