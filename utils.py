import os
import warnings
from typing import Union, List
from Bio import SeqIO
import numpy as np

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
from torch.nn.functional import log_softmax
import torch
import tempfile


snake_format = lambda s: s.replace(' ', '_').replace('-', '_').lower()


def print_if(verbose: object, thr: object, text: object) -> object:
    """
    print text if verbose > thr
    :param verbose: int
    :param thr: int
    :param text: str
    :return:
    """
    if verbose >= thr:
        print(text)


def process_fastas(text):
    """
    process multiple fasta sequences
    :return: {id: sequence}
    """
    temp = tempfile.TemporaryFile(mode='w+t')
    temp.writelines(text)
    temp.seek(0)
    ret = {seq_record.id: str(seq_record.seq) for seq_record in SeqIO.parse(temp, "fasta")}
    temp.close()
    return ret

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


def warn_if(verbose, thr, text):
    """
    print text if verbose > thr
    :param verbose: int
    :param thr: int
    :param text: str
    :return:
    """
    if verbose >= thr:
        warnings.warn(text)


def safe_post_request(session, url, timeout, verbose_level, warning_msg='connection failed', return_on_failure=None,
                      warning_thr=VERBOSE['thread_warnings'], raw_err_thr=VERBOSE['raw_warnings']):
    """
    creates a user friendly request raises warning on ConnectionError but will not crush
    verbose_level = 3 will return raw Error massage in warning
    :param session: requests session obj
    :param url: str url to query
    :param timeout: float max time to wait for response
    :param verbose_level: int
    :param warning_msg: str msg to display on failure
    :param return_on_failure: value to return upon exception
    :param raw_err_thr: int threshold to print raw error messages
    :param warning_thr: int threshold to print warning messages
    :return: response
    """
    try:
        r = session.post(url, timeout=timeout)
    except requests.exceptions.ConnectionError as e:
        warn_if(verbose_level, warning_thr, warning_msg)
        warn_if(verbose_level, raw_err_thr, f"{e}")
        return return_on_failure
    return r


def kegg_genes_in_dataset():
    return {os.path.basename(path)[:-7].replace('_', ':') for path in glob.glob(pjoin(KEGG_GENES_PATH, '*'))}


class UniprotApi:
    """
    This class is responsible to connect to online DBs and retrieve information
    """

    def __init__(self, verbose_level=1):
        self._v = verbose_level

    def fetch_uniport_sequences(self, uid):
        """
        Retrieve all known isoforms from uniprot
        :param uid: Uniprot id only primary name
        :return: {uid_iso_index: sequence}
        """
        print_if(self._v, VERBOSE['thread_progress'], "Retrieving isoforms from Uniprt...")
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        url = Q_UNIP_ALL_ISOFORMS.format(uid, uid)
        response = safe_post_request(s, url, TIMEOUT, self._v, CON_ERR_FUS.format(uid, url))
        if not response.ok:
            print_if(self._v, VERBOSE['thread_warnings'], CON_ERR_GENERAL.format('fetch_uniprot_sequences', uid))
            return {}
        if response.text == '':
            print_if(self._v, VERBOSE['thread_warnings'], f"no sequences found for {uid}")
            return {}
        return process_fastas(response.text)

    def expand_isoforms(self, ref_name, ref_mut=None, reviewed=True):
        """
        expand protein isoforms using all relevant Uniprot accession
        this will not override the default protein isoforms
        :param ref_name: protein name
        :param ref_mut: Mutation obj if given will search for isoform with the given mutation
        :return: {uid_iso_index: seq}
        """
        uids = self.uid_from_name(ref_name)['reviewed'] if reviewed else self.uid_from_name(ref_name)['all_enteries']
        isoforms = {}

        if isinstance(ref_mut, str):
            match = re.match(VARIATION_REGEX, ref_mut)
            if not match:
                raise ValueError(f"Invalid mutation format: {ref_mut}")
            orig_aa, loc, mut_aa = match.groups()
            idx = int(loc) - 1
        else:
            idx = ref_mut.loc - 1
            orig_aa = ref_mut.origAA


        # Search for isoforms that have a fitting sequence to the mutation
        for uid in uids:
            res = self.fetch_uniport_sequences(uid)
            if ref_mut:
                for iso, seq in res.items():
                    if idx >= len(seq):
                        continue
                    if seq[idx] == orig_aa:
                        return {iso: seq}
            isoforms = {**isoforms, **res}
        return isoforms if not ref_mut else {}

    def uid_from_name(self, ref_name):
        """
        return uniprot-id given a protein ref_name
        :param ref_name: protein name
        :return: {'reviewed': [...], 'non_reviewed': [...]}
        """
        ret = {'reviewed': [], 'non_reviewed': [], 'main_entery': [], 'all_enteries': [], 'aliases': []}
        query = UNIP_QUERY_URL + Q_UID_PROT_ALL.format(ref_name)
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        r = safe_get_request(s, query, TIMEOUT, self._v, CON_ERR_UFN.format(ref_name))
        if not r:
            return ret
        if r.text == '':
            return ret
        ret = self._process_uid_query(r.text)
        return ret

    @staticmethod
    def _process_uid_query(data):
        ret = {'reviewed': [], 'non_reviewed': [], 'main_entery': [], 'all_enteries': [], 'aliases': []}
        rows = data.split('\n')[1:-1]  # first is header last is blank
        if not rows:
            return ret
        main_entry = rows[0].split('\t')[UIDS_COL_IDX]  # first entry is considered main
        ret['main_entery'].append(main_entry)
        for row in rows:
            values = row.split('\t')
            entry, reviewed, gene = values[UIDS_COL_IDX], values[REVIEWED_COL_IDX], values[GENE_NAME_COL_IDX].split(" ")
            ret['all_enteries'].append(entry)
            ret['aliases'] += gene
            if reviewed == UNIP_REVIEWED:
                ret['reviewed'].append(entry)
            else:
                ret['non_reviewed'].append(entry)
        return ret


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
        :param outpath: path to save results
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
    def __init__(self, model, alphabet, testing=False):
        """
        Initialize the calculator with the model and alphabet.

        Args:
            model - the ESM model
            alphabet - the alphabet, do:
            model, alphabet = pretrained.load_model_and_alphabet(ESM1B_MODEL)
            testing - test mode: without using esm or alphabet, only statics
        """
        if not testing:
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
                repr_layers=[33],
                return_contacts=False,
            )
            logits = out['logits']
            aa_logits = logits[0, 1:-1, 4:24]  # [seq_len, 20] tensor
        return aa_logits

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
            offset = max(0, mut_idx - 511)
            sequence = sequence[offset:offset + 1022]

        batch_tokens = self._get_batch_tokens(sequence, name)
        logits = self._get_logits(batch_tokens)
        log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
        return log_probs, offset

    def score_all_mutations(self, sequence: str) -> torch.Tensor:
        """
        Scores all single amino acid substitutions for a given sequence
        using a vectorized approach after obtaining logits from the model.

        log P(mut) - log P(wt) = log-softmax(logits - wt_logit)

        Args:
            sequence: WT amino acid sequence (string)

        Returns:
            log_probs: Tensor of shape [seq_len, 20].
        """
        # --- Part 1: Sequence and Logits Acquisition (from original score_all_mutations) ---
        if len(sequence) > ESM_MAX_LENGTH:
            # using the sliding window
            logits = self.handle_long_protein(sequence)
        else:
            batch_tokens = self._get_batch_tokens(sequence, 'WT')
            logits = self._get_logits(batch_tokens)  # raw logits [seq_len, 20]

        # Ensure that the sequence length from logits matches for the WT index gathering
        effective_sequence = sequence[:logits.shape[0]]

        # --- Part 2: Vectorized Scoring ---
        # Map effective sequence to WT indices tensor
        wt_indices = torch.tensor(
            [AA_TO_INDEX_ESM.get(aa, -1) for aa in effective_sequence],
            dtype=torch.long,
            device=self.device
        )

        # Filter out invalid AAs (-1). These are positions with unknown AAs in the WT sequence.
        valid_mask = wt_indices >= 0
        valid_positions = valid_mask.nonzero(as_tuple=True)[0]

        if len(valid_positions) == 0:
            print("INFO: No valid amino acids found in the sequence for scoring.")
            return torch.empty(0, 20, device=self.device)  # Return empty tensor if no valid positions

        # For valid positions only
        wt_indices_valid = wt_indices[valid_positions]  # shape [valid_len]

        # Gather WT logits: shape [valid_len]
        wt_logits = logits[valid_positions, wt_indices_valid]

        # Subtract WT logits (broadcast): shape [valid_len, 20]
        # This is `logits - wt_logit` from the formula
        delta_logits = logits[valid_positions] - wt_logits.unsqueeze(1)

        # Apply log_softmax along AA dimension (dim=1)
        # This is `log-softmax(...)` from the formula
        log_probs = log_softmax(delta_logits, dim=1)  # shape [valid_len, 20]

        return log_probs

    def handle_long_protein(self, sequence: str) -> torch.Tensor:
        """
        Handles proteins longer than ESM_MAX_LENGTH using a sliding window.
        Each residue is assigned logits from exactly one window:
          - First window: from start to (midpoint + half stride)
          - Middle windows: exactly `stride` residues each
          - Final window (when it reaches the end): take everything until the end
        """
        seq_len = len(sequence)
        window_size = ESM_MAX_LENGTH
        stride = SLIDING_WINDOW_STRIDE
        half_window = window_size // 2
        half_stride = stride // 2

        logits_full = torch.zeros(seq_len, 20, device=self.device)

        pos = 0
        window_idx = 0

        for start in range(0, seq_len, stride):
            end = min(start + window_size, seq_len)
            subseq = sequence[start:end]  # subsequence for this window

            # Run model on this window
            batch_tokens = self._get_batch_tokens(subseq, f"WT_window_{window_idx}")
            window_logits = self._get_logits(batch_tokens)  # [window_len, 20]
            window_len = window_logits.shape[0]

            if window_idx == 0:
                # First window: take first half_window + half_stride residues
                take_len = min(half_window + half_stride, window_len)
                logits_full[0:take_len] = window_logits[0:take_len]
                pos = take_len

            elif end == seq_len:
                # Final window: take remaining residues to the end
                take_len = seq_len - pos
                logits_full[pos:] = window_logits[-take_len:]
                pos = seq_len
                break

            else:
                # Middle windows: take central `stride` residues
                center = window_len // 2
                half_stride = stride // 2
                start_idx = max(0, center - half_stride)
                end_idx = min(window_len, center + half_stride)
                take_len = end_idx - start_idx

                logits_full[pos:pos + take_len] = window_logits[start_idx:end_idx]
                pos += take_len

            window_idx += 1

        return logits_full

    @staticmethod
    def extreme_outliers_mad(data: Union[List[float], np.ndarray, pd.Series], threshold: float = 4.25) -> np.ndarray:
        """
        NaN-safe MAD outlier detection; returns indices of outliers.

        Args:
            data (Union[List[float], np.ndarray, pd.Series]): Input data with scores.
            threshold (float): The modified z-score threshold to identify an outlier.

        Returns:
            np.ndarray: An array of integer indices for the outlier values.
        """
        data_array = np.array(data, dtype=float)
        median = np.nanmedian(data_array)
        mad = np.nanmedian(np.abs(data_array - median))

        if mad == 0 or np.isnan(mad):
            return np.array([], dtype=int)

        modified_z = 0.6745 * (data_array - median) / mad
        return np.where(np.abs(modified_z) > threshold)[0]

    def normalize_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Applies the specified normalization process to the 'score' column. After
        normalization, it assigns a score of 1 to any variants that resulted in an
        empty score (e.g., mutations to a stop codon).

        Args:
            df (pd.DataFrame): The DataFrame with a 'score' column.

        Returns:
            pd.DataFrame: DataFrame with an added 'normalized_score' column.
        """
        df_processed = df.copy()
        scores = pd.to_numeric(df_processed['score'], errors='coerce')  # invalid parsing set as NaN

        # 1. Identify and mask outliers from non-NaN scores
        outlier_indices = self.extreme_outliers_mad(scores.values, threshold=4.25)
        inlier_mask = ~np.isnan(scores.values)
        if outlier_indices.size > 0:
            inlier_mask[outlier_indices] = False
        inliers = scores[inlier_mask]

        # 2. Perform min-max normalization based on inliers
        # Handle case where there are no inliers to avoid errors
        if inliers.empty:
            # If no valid scores, can't normalize; return NaN or a default value
            df_processed['normalized_score'] = np.nan
        else:
            vmin, vmax = np.nanmin(inliers), np.nanmax(inliers)

            # Avoid division by zero if all inlier values are the same
            if vmax == vmin:
                norm = pd.Series(0.5, index=scores.index)  # Assign a neutral score
            else:
                norm = (scores - vmin) / (vmax - vmin)

            # 3. Clip outliers to 0 or 1
            norm = norm.clip(lower=0.0, upper=1.0)

            # 4. Flip scores so 1 is most pathogenic
            df_processed['normalized_score'] = 1.0 - norm

        # 5. For empty scores (e.g., stop codons), assign the max pathogenicity score of 1.
        df_processed['normalized_score'] = df_processed['normalized_score'].fillna(1.0)

        return df_processed

    def save_mutation_scores_to_csv_full(self, sequence: str, input_csv_path: str, output_csv_path: str):
        """
        Computes mutation scores for a single protein and adds them to an existing SNVs CSV file.
        The scores are added to a new column named 'score'.

        Args:
            sequence: The wild-type amino acid sequence.
            input_csv_path: Path to the existing input CSV file with variant information.
            output_csv_path: Path to the CSV file where scores will be saved.
        """
        # Compute scores for the given sequence
        score_matrix = self.score_all_mutations(sequence)

        # Load the existing CSV file
        try:
            df = pd.read_csv(input_csv_path)
        except FileNotFoundError:
            raise f"Error: Input CSV file not found at {input_csv_path}"
        except Exception as e:
            raise f"Error reading CSV file: {e}"

        # Initialize a new 'score' column with NaN
        df['score'] = float('nan')

        # Iterate through the DataFrame and add scores
        for idx, row in df.iterrows():
            try:
                na_index = int(row['Start'])  # Nucleotide index
                aa_pos = na_index // 3  # Convert to codon index

                variant = str(row['Variant'])  # e.g., M0L, M2R...
                mut_aa = variant[-1]  # Extract mutant amino acid

                if mut_aa == STOP_AA:   # stop codon, so score is -inf...
                    df.at[idx, 'score'] = None      # placeholder
                    continue

                mut_idx = AA_TO_INDEX_ESM.get(mut_aa)
                if mut_idx is None or aa_pos >= score_matrix.shape[0]:    # if more sequence positions than scores calculated
                    continue

                df.at[idx, 'score'] = score_matrix[aa_pos, mut_idx].item()
            except Exception as e:
                print(f"[Warning] Could not assign score for row {idx} in {input_csv_path}: {e}")

        # Normalize scores by removing outliers and using min-max normalization
        df_normalized = self.normalize_scores(df)

        # Save the updated DataFrame to a new CSV file
        df_normalized.to_csv(output_csv_path, index=False)

    def save_mutation_scores_to_csv(self, sequence: str, csv_path: str):
        """In place alias of save_mutation_scores_to_csv_full."""
        return self.save_mutation_scores_to_csv_full(sequence, csv_path, csv_path)

