from definitions import *
from bravado.client import SwaggerClient

def set_cbio_api_call():
    cbioportal = SwaggerClient.from_url(CBIO_API_URL,
                                        config={"validate_requests": False, "validate_responses": False,
                                                "validate_swagger_spec": False})
    for a in dir(cbioportal):
        cbioportal.__setattr__(a.replace(' ', '_').lower(), cbioportal.__getattr__(a))
    return cbioportal

def download_study_mutations(cbio_api, study):
    muts = cbio_api.mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
        molecularProfileId=f"{study}_mutations",
        # {study_id}_mutations gives default mutations profile for study
        sampleListId=f"{study}_all",  # {study_id}_all includes all samples
        projection="DETAILED"  # include gene info
    return muts

