import time
from utils import *
from Kegg import *

def init_kegg_genome(recalc=False):
    """Initialize KEGG genome by creating KeggGene objects for all genes in KEGG database."""
    genes = KeggApi().get_all_genes().keys()
    if not recalc:
        genes = set(genes) - kegg_genes_in_dataset()
    genes = list(genes)
    def objects_creator(*gene_ids):
        kegg = KeggApi()
        genes_dict = kegg.genes_info(list(gene_ids))
        for kegg_id, data in genes_dict.items():
            gene = KeggGene(kegg_id, default_init=False)
            gene.create_from_dict(data)

    target = objects_creator
    #  at most 10 genes per request
    tasks = [genes[i: i + KEGG_MAX_IDS] for i in range(0, len(genes), KEGG_MAX_IDS)]
    multiprocess_task(tasks=tasks, target=target, workers=KEG_API_RECOMMENDED_WORKERS)


if __name__ == '__main__':
    # Step 1: Initialize KEGG genome (download all genes)
    init_kegg_genome()

    # Step 2: Get all pathways and modules
    kegg = KeggApi()
    pathways = kegg.get_all_pathways()
    modules = kegg.get_all_modules()

    # Step 3: Generate SNV CSVs for each pathway
    for pathway_id in pathways:
        print(f"Processing pathway: {pathway_id}")
        start = time.time()
        network = KeggNetwork(pathway_id, 'pathway')
        network.all_snvs()  # will save to default path
        print(f"Finished processing {pathway_id} in {time.time() - start:.2f} seconds")

    # Step 4: Generate SNV CSVs for each module
    for module_id in modules:
        print(f"Processing module: {module_id}")
        network = KeggNetwork(module_id, 'module')
        network.all_snvs()  # will save to default path

