import os
import warnings
from definitions import *
from bravado.client import SwaggerClient
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool as Pool
import copy
import re
import glob
import torch.nn.functional as tnf
import torch
from typing import Dict
import csv


snake_format = lambda s: s.replace(' ', '_').replace('-', '_').lower()


def read_in_chunks(array, chunk_size):
    for i in range(0, len(array), chunk_size):
        yield array[i:i + chunk_size]


def save_dict(data, path):
    with open(path, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def load_dict(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


def save_obj(obj, path):
    with open(path, "wb") as file:
        pickle.dump(obj.__dict__, file)


def load_obj(obj, path, name=''):
    with open(path, 'rb') as file:
        if os.path.getsize(path) > 0:
            obj.__dict__ = pickle.load(file)
        else:
            raise NameError(LOAD_OBJ_ERROR.format(name))


def create_session(header, retries=5, wait_time=0.5, status_forcelist=None):
    """
    Creates a session using pagination
    :param header: str url header session eill apply to
    :param retries: int number of retries on failure
    :param wait_time: float time (sec) between attempts
    :param status_forcelist: list HTTP status codes that we should force a retry on
    :return: requests session
    """
    s = requests.Session()
    retries = Retry(total=retries,
                    backoff_factor=wait_time,
                    status_forcelist=status_forcelist)

    s.mount(header, HTTPAdapter(max_retries=retries))
    return s


def safe_get_request(session, url, timeout=TIMEOUT, warning_msg='connection failed', return_on_failure=None):
    """
    creates a user friendly request raises warning on ConnectionError but will not crush
    verbose_level = 3 will return raw Error massage in warning
    :param session: requests session obj
    :param url: str url to query
    :param timeout: float max time to wait for response
    :param warning_msg: str msg to display on failure
    :param return_on_failure: value to return upon exception
    :return: response
    """
    try:
        r = session.get(url, timeout=timeout)
    except requests.exceptions.ConnectionError as e:
        warnings.warn(warning_msg)
        return return_on_failure
    return r


def multiprocess_task(tasks, target, workers=None, callback=lambda x: x):
    """
    :param tasks: list of iterables
    :param target: callable
    :param workers: optional int number of CPUs that will be used otherwise maximum available
    :param callback: callable
    :return:
    """
    workers = workers if workers else cpu_count()
    with Pool(workers) as p:
        try:
            for status in p.starmap(target, tasks):
                callback(status)
        except Exception as e:
            print(f"Multiprocessing task failed: {e}")


def kegg_genes_in_dataset():
    return {os.path.basename(path)[:-7].replace('_', ':') for path in glob.glob(pjoin(KEGG_GENES_PATH, '*'))}

class CbioApi:
    """api for cbio portal"""
    # cbio is portal for cancer genomics data

    def __init__(self, ):
        """Constructor for Cbio"""
        cbioportal = SwaggerClient.from_url(CBIO_API_URL,
                                            config={"validate_requests": False, "validate_responses": False,
                                                    "validate_swagger_spec": False})
        for a in dir(cbioportal):
            cbioportal.__setattr__(snake_format(a), cbioportal.__getattr__(a))
        self.api = cbioportal

    @staticmethod
    def set_cbio_api_call():
        cbioportal = SwaggerClient.from_url(CBIO_API_URL,
                                            config={"validate_requests": False, "validate_responses": False,
                                                    "validate_swagger_spec": False})
        for a in dir(cbioportal):
            cbioportal.__setattr__(snake_format(a), cbioportal.__getattr__(a))
        return cbioportal

    def get_sample_cancer_type(self, study_id, sample_id):
        data = self.api.Clinical_Data.getAllClinicalDataOfSampleInStudyUsingGET(studyId=study_id,
                                                                                sampleId=sample_id,
                                                                                attributeId='CANCER_TYPE')
        return snake_format(data.result()[0].value)

    def download_study_mutations(self, study):
        muts = self.api.mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
            molecularProfileId=f"{study}_mutations",
            # {study_id}_mutations gives default mutations profile for study
            sampleListId=f"{study}_all",  # {study_id}_all includes all samples
            projection="DETAILED")  # include gene info
        return muts

    @staticmethod
    def study_to_csv(results, outpath='', remove_duplicates=True):
        """
        :param remove_duplicates:
        :param outpath:
        :param results: bravado.http_future.HttpFuture object
        :return: csv in FamAnalysis format
        """
        mutations = results.result()
        data = [(m.chr, m.startPosition, m.endPosition, m.referenceAllele, m.variantAllele, m.gene.hugoGeneSymbol,
                 m.proteinChange, m.patientId, m.uniquePatientKey, m.sampleId, m.studyId, m.ncbiBuild)
                for m in mutations if m.mutationType == MISSENSE_MUTATION]
        df = pd.DataFrame.from_records(data, columns=STUDY_COLUMNS)
        #  drop duplicate mutations of the same patient
        #  mutations can repeat in the same patient in the same study if there are multiple samples per patient
        if remove_duplicates:
            df.drop_duplicates(keep='first', inplace=True, ignore_index=True, subset=DUPLICATE_EXCLUSION_COLUMNS)
        if outpath:
            df.to_csv(outpath)
        return df

    def cancer_types_dict(self):
        """
        :return: dict {cancer_type : cbio_short_name}
        """
        all_types = self.api.Cancer_Types.getAllCancerTypesUsingGET().result()
        return {snake_format(cancer_type.name): cancer_type.shortName for cancer_type in all_types}

    def all_studies_by_keyword(self, keyword, outpath=''):
        """
        :param keyword: abbreviated cancer type
        :return: all studies with samples of cancer_type == keyword
        """
        studies = self.api.Studies.getAllStudiesUsingGET(keyword=keyword).result()
        #  make sure cancer type is correct
        studies = [study for study in studies if study.cancerTypeId == keyword]
        study_ids = [study.studyId for study in studies]
        study_names = [study.name for study in studies]
        if outpath:
            with open(outpath, 'w+') as file:
                for id, name in zip(study_ids, study_names):
                    file.write(f"{name} \t {id}\n")
        return study_ids, study_names


class KeggApi:
    """api for kegg"""

    def __init__(self):
        """
        constructor for Kegg
        """
        self.api = create_session(KEGG_API_URL, retries=RETRIES, wait_time=WAIT_TIME,
                                  status_forcelist=RETRY_STATUS_LIST)

    @staticmethod
    def format_multiple_genes(genes):
        return '+'.join(genes)

    def _kegg_command(self, command, *params, verbose=True):
        """
        Kegg general command
        :param database: pathway | brite | module | ko | <org> | vg | vp | ag | genome | compound |
             glycan | reaction | rclass | enzyme | network | variant | disease |
             drug | dgroup | disease_ja | drug_ja | dgroup_ja | compound_ja
        :param species: str
        :param verbose: bool if true request will be printed
        :return:
        """
        query = command.format(*params)
        if verbose:
            print(command.format(*params))
        data = safe_get_request(self.api, query)
        if not data:
            raise ConnectionError(f'querry: {query} --> Failed')
        if not data.ok:
            raise ConnectionError(f'querry: {query} --> Failed with code {data.status_code}')
        return data.text

    def kegg_command(self, command_type, *params):
        """
        Factory for kegg commands
        :param command_type: one of list | link | conv | get
        :param params: parameters for query
        :return:
        """
        assert command_type in COMMAND_TYPES, 'command_type must be one of link | list | conv | get'
        if command_type == 'list':
            return self._kegg_command(KEGG_LIST_COMMAND, *params)
        elif command_type == 'link':
            return self._kegg_command(KEGG_LINK_COMMAND, *params)
        elif command_type == 'conv':
            return self._kegg_command(KEGG_CONV_COMMAND, *params)
        elif command_type == 'get':
            return self._kegg_command(KEGG_GET_COMMAND, *params)

    def _chunk_request(self, kegg_ids, kegg_command, *params):
        """
        efficiently process request of multiple entries into chunks
        :param kegg_command: kegg command type one of list | link | conv | get
        :param kegg_ids: str or list of kegg ids
        :param params: optional extra parameters for command
        :return: raw data of all kegg_ids with the given the kegg_command
        """
        all_data = ''
        if isinstance(kegg_ids, str):
            kegg_ids = [kegg_ids]
        #  kegg takes 10 values at a time
        chunks = [kegg_ids[i:i + KEGG_MAX_IDS] for i in range(0, len(kegg_ids), KEGG_MAX_IDS)]
        for chunk in chunks:
            query = self.format_multiple_genes(chunk)
            data = self.kegg_command(kegg_command, query, *params)
            all_data += data
        return all_data

    def get_all_genes(self, species=KEGG_HOMO_SAPIENS):
        """
        :param species: str
        :return: list of Kegg <species> genes
        """
        data = self.kegg_command('list', species, '')
        return self._process_response(data)

    def get_all_pathways(self, species=KEGG_HOMO_SAPIENS):
        """
        :param species: str
        :return: list of Kegg <species> pathways
        """
        data = self.kegg_command('list', 'pathway', species)
        return self._process_response(data)

    def get_all_modules(self):
        """
        :return: list of Kegg <species> pathways
        """
        data = self.kegg_command('list', 'module', '')
        return self._process_response(data)

    def module_orthologs(self, module_id):
        """
        list of ortholog involved in kegg module
        :param module_id:
        :return:
        """
        data = self.kegg_command('link', 'ko', module_id)
        orthologs = {ortholog.split(':')[1] for ortholog in self._process_response(data, return_as_set=True)}
        return orthologs

    def ortholog_genes(self, orthologs, species=KEGG_HOMO_SAPIENS):
        """
        :type orthologs: set or string
        :param species:
        :return:
        """
        if isinstance(orthologs, (set, list)):
            orthologs = self.format_multiple_genes(orthologs)
        data = self.kegg_command('link', species, orthologs)
        return self._process_response(data, return_as_set=True)

    @staticmethod
    def _process_response(pathway_data, return_as_list=False, return_as_set=False):
        if return_as_set:  # correction to support older version
            return_as_list = True
        res = {} if not return_as_list else []
        rows = pathway_data.split('\n')
        for row in rows[:-1]:
            items = row.split('\t')
            kegg_id = items[0]
            desc = items[1] if len(items) >= 2 else ''
            if return_as_list:
                res.append(desc)
            else:
                res[kegg_id] = desc
        if return_as_set:
            res = set(res)
        return res

    def get_pathway_info(self, pathway_id):
        """
        retrieves raw pathway information in xml format
        :param pathway_id: str - kegg pathway id
        :return:
        """
        return self.kegg_command('get', pathway_id, 'kgml')

    def get_pathway_gene_list(self, pathway_id):
        """
        retrieves set of genes in a pathway
        :param pathway_id: str kegg pathway id
        :return:
        """
        data = self.kegg_command('link', KEGG_HOMO_SAPIENS, pathway_id)
        genes = self._process_response(data, return_as_set=True)
        return genes

    def get_module_gene_list(self, module_id):
        """
        retrieves set of genes in a module
        :param module_id: str kegg pathway id
        :return:
        """
        orthologs = self.module_orthologs(module_id)
        return self.ortholog_genes(orthologs)

    def get_gene_list(self, kegg_id):
        """
        given general kegg_id retrieves all genes in network
        :param kegg_id:
        :return: set list of kegg gene ids
        """
        if kegg_id.startswith(KEGG_PATHWAY_PREFIX):
            return self.get_pathway_gene_list(kegg_id)
        elif kegg_id.startswith(KEGG_MODULE_PREFIX):
            return self.get_module_gene_list(kegg_id)
        else:
            raise ValueError(NETWORK_ID_ERROR)

    def convert_gene_names(self, genes, database='uniprot'):
        """
        convert gene names between datasets
        :param genes: str or list of gene names
        :param database: one of genes | ncbi-geneid | ncbi-proteinid | uniprot
        :return: dict containing primary entry and secondary entries
        """
        res = {'primary': None, 'secondary': set()}
        if isinstance(genes, (set, list)):
            genes = self.format_multiple_genes(genes)
        data = self.kegg_command('conv', database, genes)
        genes = self._process_response(data, return_as_set=False, return_as_list=True)
        if (not genes) or (genes == EMPTY_LIST):
            return res
        #  first instance is saved as primary
        res['primary'] = genes[0].split(':')[1]
        for gene in genes[1:]:
            res['secondary'].add(gene.split(':')[1])
        return res

    def gene_seq(self, genes, seq_type='aaseq'):
        """
        retrieves amino-acid or dna sequence of gene
        :param genes: str or list of gene names
        :param seq_type: one of aaseq | ntseq
        :return: dict {kegg_id : seq}
        """
        data = self._chunk_request(genes, 'get', seq_type)
        return self._gene_seq_to_dict(data)

    def genes_info(self, genes):
        """
        :param genes: str or list of gene names
        :return: dict {gene_id : dict gene details}
        Note: in rare cases some genes will be skipped
        """
        data = self._chunk_request(genes, 'get', '')
        return self._process_genes_info(data)

    @staticmethod
    def _gene_seq_to_dict(data):
        res = {}
        rows = data.split('\n')[:-1]
        gene_name = ''
        for row in rows:
            if row.startswith('>'):  # header
                gene_name = row.split(' ')[0][1:]
                if gene_name not in res:
                    res[gene_name] = ''
                continue
            else:
                res[gene_name] += row
        return res

    def _process_genes_info(self, data):
        genes_dict = {}
        genes = data.split(GENE_SEPERATOR)[:-1]
        for gene in genes:
            data = self._process_single_gene(gene)
            if not data:
                continue
            genes_dict = genes_dict | data
        return genes_dict

    @staticmethod
    def _process_single_gene(data):
        # define default gene data
        gene_data = copy.deepcopy(GENE_DATA)
        kegg_id = None
        aa_seq_flag, na_seq_flag = False, False
        aa_seq_len, na_seq_len = None, None

        # process each row in the gene data
        for row in data.split('\n'):
            if not row.strip():
                continue
            values = row.split()
            # while flags are True len(values) must be 1
            if len(values) > 1:
                aa_seq_flag, na_seq_flag = False, False
            # append sequences until flag is False
            if aa_seq_flag:
                gene_data['aa_seq'] += values[0]
                continue
            if na_seq_flag:
                gene_data['na_seq'] += values[0]
                continue
            try:
                title = values[0]
            except IndexError:  # happens in rare cases [empty lines, corrupted data...]
                print(f"[KEGG PARSE ERROR] Malformed line in gene block: {row}")
                continue
            # process each title
            if title == 'ENTRY':
                kegg_id = f'hsa:{values[1]}'
                gene_data['coding_type'] = values[2]
            if title == 'SYMBOL':
                gene_data['ref_names'] = values[1:]
            if title == 'POSITION':
                line_data = values[1].split(':')
                gene_data['chr'] = line_data[0]
                try:
                    res = re.search(KEG_POSITION_RE, line_data[1])
                    gene_data['start'], gene_data['end'] = int(res.groups()[0]), int(res.groups()[1])
                except IndexError:
                    gene_data['chr'] = values[1]
            if title.startswith('UniProt'):
                ids = values[1:]
                gene_data['uniprot_id'] = { 'primary': ids[0] if ids else None,
                                            'secondary': set(ids[1:]) if len(ids) > 1 else set()  }
            if title == 'AASEQ':
                aa_seq_flag = True
                aa_seq_len = int(values[1])
                continue
            if title == 'NTSEQ':
                na_seq_flag = True
                na_seq_len = int(values[1])
                continue
        # check if all required fields are present
        if gene_data['aa_seq']:
            assert len(gene_data['aa_seq']) == aa_seq_len, 'aa_seq length does not match ' + kegg_id
        if gene_data['na_seq']:
            assert len(gene_data['na_seq']) == na_seq_len, 'na_seq length does not match ' + kegg_id
        assert kegg_id, 'Unable to find Kegg id'
        gene_data['kegg_id'] = kegg_id
        return {kegg_id: gene_data}



class WildtypeMarginalsCalculator:
    def __init__(self, model, alphabet):
        """
        Initialize the calculator with the model and alphabet.
        """
        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.model = model.to(self.device)
        self.alphabet = alphabet
        self.tokenizer = alphabet.get_batch_converter()

    def _get_batch_tokens(self, sequence: str, name='WT') -> torch.Tensor:
        """
        Tokenizes the input sequence using the ESM tokenizer.
        """
        batch = [(name, sequence)]
        _, _, tokens = self.tokenizer(batch)
        return tokens.to(self.device)

    def _get_logits(self, batch_tokens: torch.Tensor) -> torch.Tensor:
        """
        Returns raw logits from the model for amino acid predictions.
        """
        self.model.eval()
        with torch.no_grad():
            out = self.model(
                batch_tokens,
                repr_layers=REP_LAYERS,
                return_contacts=False,
            )
            logits = out['logits']
            # remove [CLS] and [EOS] tokens
            if ESM1B_MODEL == 'esm1_t6_43M_UR50S':
                aa_logits = logits[0, 1:, 4:24]
            else:
                aa_logits = logits[0, 1:-1, 4:24]  # [seq_len, 20]
        return aa_logits

    @staticmethod
    def esm_process_long_sequences(seq, loc):
        """
        trims seq to len < 1024 using mut Object. to fit for bert model
        :return: offset from orig location, trimmed sequence
        """
        thr = ESM_MAX_LENGTH // 2
        left_bound = 0 if loc - thr < 0 else loc - thr
        right_bound = len(seq) if loc + thr > len(seq) else loc + thr
        left_excess = 0 if loc - thr > 0 else abs(loc - thr)
        right_excess = 0 if loc + thr <= len(seq) else loc + thr - len(seq)
        if (left_excess == 0) and (right_excess == 0):
            return left_bound, seq[left_bound:right_bound]
        if left_excess > 0:
            return 0, seq[left_bound:right_bound + left_excess]
        if right_excess > 0:
            return left_bound - right_excess, seq[left_bound - right_excess:right_bound]

    def compute_log_marginals(self, sequence: str, mut_idx: int = None, name='WT'):
        """
        Compute log-softmax over amino acid logits to get log-marginals.
        Handles long sequences by slicing around mutation index.

        Returns:
            log_probs: Tensor of shape [seq_len, 20]
            offset: offset applied to sequence (for long seqs); 0 if no slicing
        """
        offset = 0
        if len(sequence) > ESM_MAX_LENGTH:
            assert mut_idx is not None, "Mutation index required for long sequence processing."
            offset, sequence = self.esm_process_long_sequences(sequence, mut_idx)

        batch_tokens = self._get_batch_tokens(sequence, name)
        logits = self._get_logits(batch_tokens)
        log_probs = tnf.log_softmax(logits, dim=-1)
        return log_probs, offset

    @staticmethod
    def score_mutation(log_probs: torch.Tensor, aa_mut, offset: int = 0) -> float:
        """
        Compute mutation score from log-marginals (log p(mut) - log p(wt)).

        Args:
            log_probs: Tensor [seq_len, 20]
            aa_mut: Mutation object
            offset: offset applied if sequence was sliced

        Returns:
            float: mutation log-ratio score
        """
        idx = aa_mut.mut_idx - offset
        wt_idx = AA_TO_INDEX_ESM[aa_mut.wt_aa]
        mut_idx = AA_TO_INDEX_ESM[aa_mut.change_aa]
        return (log_probs[idx, mut_idx] - log_probs[idx, wt_idx]).item()

    @staticmethod
    def score_all_mutations(log_probs: torch.Tensor, wt_sequence: str, offset: int = 0) -> Dict[str, float]:
        """
        Score all possible single mutations at all positions in the sequence using vectorized operations.

        Args:
            log_probs: Tensor of shape [seq_len, 20], log-marginals from WT
            wt_sequence: Full WT sequence (untruncated)
            offset: Offset from slicing, if applicable

        Returns:
            Dict[str, float]: Mutation scores in format { "R29L": score, ... }
        """
        seq_len = log_probs.size(0)
        device = log_probs.device
        mutation_scores = {}

        # Map WT AAs in current window to indices
        wt_aas = [wt_sequence[i + offset] for i in range(seq_len)]
        valid_mask = torch.tensor([aa in AA_TO_INDEX_ESM for aa in wt_aas], dtype=torch.bool, device=device)
        wt_indices = torch.tensor(
            [AA_TO_INDEX_ESM[aa] if aa in AA_TO_INDEX_ESM else -1 for aa in wt_aas],
            dtype=torch.long,
            device=device
        )

        # Prepare tensor for all mutant AA indices
        all_mutant_aas = [aa for aa in VALID_AA]
        all_mutant_indices = torch.tensor(
            [AA_TO_INDEX_ESM[aa] for aa in all_mutant_aas],
            dtype=torch.long,
            device=device
        )  # [20]

        for i in range(seq_len):
            if not valid_mask[i]:
                continue
            wt_aa = wt_aas[i]
            wt_idx = wt_indices[i].item()
            wt_log_probs = log_probs[i]  # shape: [20]

            # Compute delta = log p(mut_aa) - log p(wt_aa)
            delta = log_probs[i][all_mutant_indices] - wt_log_probs[wt_idx]  # [20]

            for j, mut_aa in enumerate(all_mutant_aas):
                if mut_aa == wt_aa:
                    continue  # skip WT->WT
                mutation_scores[f"{wt_aa}{i + offset + 1}{mut_aa}"] = delta[j].item()

        return mutation_scores

    @staticmethod
    def save_mutation_scores_to_csv(
            mutation_scores: Dict[str, float],
            output_path: str,
            protein_id: str,
    ):
        """
        Save mutation scores to a CSV file in the format:
        ,Chr,Start,End,Ref,Alt,Protein,Variant

        Args:
            mutation_scores: Dictionary like { "M1L": score, ... }
            output_path: Path to the CSV file
            protein_id: e.g. 'P29401'
        """
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Protein', 'Variant'])

            for idx, (variant, score) in enumerate(mutation_scores.items()):
                wt_aa = variant[0]
                pos = int(variant[1:-1])
                mut_aa = variant[-1]
                writer.writerow([
                    idx,
                    '-',  # Chr placeholder
                    pos - 1,  # Start (0-based)
                    pos - 1,  # End
                    wt_aa.lower(),  # Ref nucleotide placeholder (dummy)
                    mut_aa.lower(),  # Alt
                    protein_id,
                    variant
                ])

    def save_all_proteins_to_csv(
            self,
            protein_sequences: Dict[str, str],
            output_path: str
    ):
        """
        Compute mutation scores for multiple proteins and save all results to one CSV file.

        Args:
            protein_sequences: Dict mapping protein IDs to sequences.
            output_path: Path to save the combined CSV file. [SCORING_MUTATIONS_PATH]
        """
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Protein', 'Variant'])

            row_idx = 0
            for protein_id, sequence in protein_sequences.items():
                log_probs, offset = self.compute_log_marginals(sequence)
                mutation_scores = self.score_all_mutations(log_probs, sequence, offset)

                for variant, score in mutation_scores.items():
                    wt_aa = variant[0]
                    pos = int(variant[1:-1])
                    mut_aa = variant[-1]
                    writer.writerow([
                        row_idx,
                        '-',  # Chr placeholder
                        pos - 1,  # Start (0-based)
                        pos - 1,  # End
                        wt_aa.lower(),  # Ref (dummy nucleotide)
                        mut_aa.lower(),  # Alt (dummy nucleotide)
                        protein_id,
                        variant
                    ])
                    row_idx += 1


