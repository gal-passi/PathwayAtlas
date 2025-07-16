import pickle
from os.path import join as pjoin
import re


def save_dict(data, path):
    with open(path, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def load_dict(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


#  DIRECTORIES

DB = 'data'
CBIO_PATH = pjoin(DB, 'cbio')
KEGG_PATH = pjoin(DB, 'kegg')
STUDIES_PATH = pjoin(CBIO_PATH, 'studies')
KEGG_GENES_PATH = pjoin(KEGG_PATH, 'genes')
KEGG_PATHWAYS_PATH = pjoin(KEGG_PATH, 'pathways')
KEGG_PATHWAY_OBJECTS_PATH = pjoin(KEGG_PATHWAYS_PATH, 'objects')
KEGG_PATHWAY_MUTATIONS_PATH = pjoin(KEGG_PATHWAYS_PATH, 'snvs')

DIRS_TO_CREATE = [DB, CBIO_PATH, KEGG_PATH, STUDIES_PATH, KEGG_GENES_PATH, KEGG_PATHWAYS_PATH,
                  KEGG_PATHWAY_OBJECTS_PATH,KEGG_PATHWAY_MUTATIONS_PATH]

#   REQUESTS AND OS CONSTANTS

TIMEOUT = 3.0
WAIT_TIME = 1.0
RETRIES = 7
RETRY_STATUS_LIST = [429, 500, 502, 503, 504, 403, 400]
DEFAULT_HEADER = "https://"
WORKERS = None  # use default amount of CPUs
KEG_API_RECOMMENDED_WORKERS = 6

#  BIOLOGIC CONSTANTS

NA_COUPLE = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
NA_CHANGE = {'a': 'tcg', 't': 'acg', 'c': 'gta', 'g': 'cta'}
STOP_CODONS = ['tag', 'taa', 'tga']
STOP_AA = '_'
CODON_LENGTH = 3
CODON_TRANSLATOR = {'ata': 'I', 'atc': 'I', 'att': 'I', 'atg': 'M', 'aca': 'T',
                    'acc': 'T', 'acg': 'T', 'act': 'T', 'aac': 'N', 'aat': 'N',
                    'aaa': 'K', 'aag': 'K', 'agc': 'S', 'agt': 'S', 'aga': 'R',
                    'agg': 'R', 'cta': 'L', 'ctc': 'L', 'ctg': 'L', 'ctt': 'L',
                    'cca': 'P', 'ccc': 'P', 'ccg': 'P', 'cct': 'P', 'cac': 'H',
                    'cat': 'H', 'caa': 'Q', 'cag': 'Q', 'cga': 'R', 'cgc': 'R',
                    'cgg': 'R', 'cgt': 'R', 'gta': 'V', 'gtc': 'V', 'gtg': 'V',
                    'gtt': 'V', 'gca': 'A', 'gcc': 'A', 'gcg': 'A', 'gct': 'A',
                    'gac': 'D', 'gat': 'D', 'gaa': 'E', 'gag': 'E', 'gga': 'G',
                    'ggc': 'G', 'ggg': 'G', 'ggt': 'G', 'tca': 'S', 'tcc': 'S',
                    'tcg': 'S', 'tct': 'S', 'ttc': 'F', 'ttt': 'F', 'tta': 'L',
                    'ttg': 'L', 'tac': 'Y', 'tat': 'Y', 'taa': '_', 'tag': '_',
                    'tgc': 'C', 'tgt': 'C', 'tga': '_', 'tgg': 'W'}

#  FAMANALYSIS

FAMANALYSIS_COLUMNS = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Protein', 'Variant']

#   CBIOPORTAL

CBIO_API_URL = 'https://www.cbioportal.org/api/v2/api-docs'
MISSENSE_MUTATION = 'Missense_Mutation'
CANCER_TYPES_DIR = pjoin(CBIO_PATH, 'cbio_cancer_types.pickle')
CBIO_CANCER_TYPES = load_dict(CANCER_TYPES_DIR)
STUDY_COLUMNS = FAMANALYSIS_COLUMNS + ['PatientId', 'PatientKey', 'SampleId', 'StudyId', 'RefDNA']
# exclude only on patient key and protein change to avoid problems with hg19/hg18
DUPLICATE_EXCLUSION_COLUMNS = FAMANALYSIS_COLUMNS + ['PatientKey']

#  KEGG
KEGG_HSA_PATHWAYS_DIR = pjoin(KEGG_PATH, 'kegg_hsa_pathways.pickle')  # {pathway_id : desc}

KEGG_API_URL = 'https://rest.kegg.jp'
COMMAND_TYPES = ['link', 'list', 'conv', 'get']
KEGG_LIST_COMMAND = KEGG_API_URL + '/list/{}/{}'
KEGG_LINK_COMMAND = KEGG_API_URL + '/link/{}/{}'
KEGG_CONV_COMMAND = KEGG_API_URL + '/conv/{}/{}'
KEGG_GET_COMMAND = KEGG_API_URL + '/get/' + '{}/{}'

KEGG_GENE_SEQ_URL = f'{KEGG_API_URL}/get/' + '{}/{}'
KEGG_HOMO_SAPIENS = 'hsa'
KEGG_PATHWAY_PREFIX = 'hsa'
KEGG_MODULE_PREFIX = 'M'

KEGG_MAX_IDS = 10
NETWORK_TYPES = ('module', 'pathway')
EMPTY_SET = {''}
EMPTY_LIST = ['']
GENE_SEPERATOR = '///\n'

GENE_DATA = {'kegg_id': None, 'uniprot_id': None, 'aa_seq': '', 'na_seq': '',
             'chr': None, 'start': None, 'end': None, 'coding_type': None, 'ref_names': None}

KEG_POSITION_RE = "(\d+)\.{2}(\d+)"

#  ERRORS

NETWORK_TYPE_ERROR = f'Network type must be one of: {", ".join(NETWORK_TYPES)}'
NETWORK_ID_ERROR = f'KEGG id must be of a KEGG module or KEGG pathway'
LOAD_OBJ_ERROR = 'Data missing or invalid for {}. ' \
                 '\nDelete instance from DB and recreate the object'
