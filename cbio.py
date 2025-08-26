from utils import *
from Kegg import *
import os

"""
This module provides functions to interact with the cBioPortal API.
Run the `get_studies` function to download mutation data for all studies.
"""


def check_for_duplicates(dfs):
    """Check for duplicate mutations in all the studies.
        @param dfs: dictionary of study id to df of that study, returned by `get_studies()`.
    """
    all_patients = []
    for study, df in dfs.items():
        patients = df['PatientId'].unique()
        all_patients.extend(patients)

    all_df = pd.DataFrame(all_patients, columns=['PatientId', 'StudyId'])

    overlaps = (all_df.groupby('PatientId')['StudyId']
                .nunique()
                .reset_index()
                .query("StudyId > 1")
                )

    if overlaps.empty:
        print("No overlapping patients found.")
    else:
        print(f"found {len(overlaps)} overlapping patients.")


def get_studies(cbio: CbioApi):
    """Retrieve all cBioPortal studies.
        Downloads mutations for each study and saves them as CSV files.
        Returns a dictionary of study_id to DataFrame of mutations.
    """
    dfs = {}
    all_studies = cbio.api.Studies.getAllStudiesUsingGET().result()

    for study in all_studies:
        study_id = study.studyId
        try:
            # Download mutations
            results = cbio.download_study_mutations(study_id)

            # Define file name for each study
            outpath = f"data/cbio/studies/{study_id}_mutations.csv"

            if os.path.exists(outpath):
                print(f"{study_id} mutations already downloaded.")
                dfs[study_id] = pd.read_csv(outpath)
            else:
                # Save DataFrame to CSV
                df = cbio.study_to_csv(results, outpath=outpath)
                dfs[study_id] = df
                print(f"Saved {study_id}_mutations.csv")
        except Exception as e:
            print(f"Skipping {study_id}: {e}")

    print("Downloaded all studies.")
    return dfs


def merge_studies(cbio, dfs, remove_duplicates=True):
    """
    merge the studies of the same cancer type into a single dataframe.
    @param cbio: cBioPortal API object
    @param dfs: dictionary of study id to DataFrame of that study, returned by 'get_studies()'.
    @param remove_duplicates: remove duplicate studies if True.
    """
    cancer_dfs = {}
    cancer_types_dict = cbio.cancer_types_dict()

    for cancer_type, short_cancer_name in cancer_types_dict.items():

        # Get all studies for the given cancer type
        study_ids, study_names = cbio.all_studies_by_keyword(short_cancer_name.lower())
        study_ids = [id for id in study_ids if id in dfs.keys()]

        if not study_ids:
            print(f"No studies found for cancer type {short_cancer_name.lower()}")
            continue  # skip to the next cancer type

        # Merge them all to one Dataframe
        merged_df = pd.concat([dfs[study_id] for study_id in study_ids], ignore_index=True)

        # Optionally remove duplicates across studies
        if remove_duplicates:
            merged_df.drop_duplicates(keep='first', inplace=True, ignore_index=True, subset=DUPLICATE_EXCLUSION_COLUMNS)

        # Define file name for each cancer
        outpath = f"data/cbio/cancers/{short_cancer_name.lower()}_mutations.csv"

        if os.path.exists(outpath):
            print(f"{cancer_type} mutations already downloaded.")
            cancer_dfs[cancer_type] = pd.read_csv(outpath)
        else:
            # Save merged DataFrame to CSV
            merged_df.to_csv(outpath, index=False)
            cancer_dfs[cancer_type] = merged_df
            print(f"Saved {short_cancer_name.lower()}_mutations.csv")

    print("Merged all cancer studies.")
    return cancer_dfs


def get_uniprot_canonical_sequence(gene_name: str, organism: str = "Homo sapiens") -> str:
    """
    Fetches the canonical UniProt protein sequence for a given gene name and species.

    Args:
        gene_name (str): Gene symbol (e.g., "DNMT1").
        organism (str): Organism name (default = "Homo sapiens").

    Returns:
        str: The canonical protein sequence (FASTA format).
    """
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f'gene_exact:{gene_name} AND organism_name:"{organism}" AND reviewed:true',
        "fields": "accession",
        "size": 1
    }
    response = requests.get(search_url, params=params)
    response.raise_for_status()
    data = response.json()

    if not data.get("results"):
        raise ValueError(f"No reviewed UniProt entry found for gene {gene_name} in {organism}")

    accession = data["results"][0]["primaryAccession"]

    # Step 2: Fetch the FASTA sequence for that accession
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    fasta_response = requests.get(fasta_url)
    fasta_response.raise_for_status()
    fasta_text = fasta_response.text

    # Step 3: Parse FASTA to get sequence string
    lines = fasta_text.splitlines()
    sequence = "".join(lines[1:])  # skip header
    return sequence


def check_sequences(cancer_dfs):

    for cancer, df in cancer_dfs.items():
        for row in df.itertuples(index=False):
            protein = row.Protein
            variant = row.Variant

            sequence = get_uniprot_canonical_sequence(protein)
            match = re.match(r"([A-Z])(\d+)([A-Z])", variant)
            if not match:
                raise ValueError(f"Invalid mutation format: {variant}")

            ref_aa, pos, alt_aa = match.groups()
            pos = int(pos)

            # UniProt is 1-based indexing, Python is 0-based
            # TODO: check if this works now
            if pos > len(sequence) or pos < 1:
                print(
                    f"Variant position {pos} is out of range for protein {protein} (length {len(sequence)}). Skipping.")
                continue

            seq_aa = sequence[pos - 1]
            if seq_aa != ref_aa:
                print(f"Wrong sequence for protein {protein}: expected {ref_aa} at position {pos}, found {seq_aa}")

        print(f"Finished validating sequences in cancer type {cancer}.")
    print("Finished validating sequences.")