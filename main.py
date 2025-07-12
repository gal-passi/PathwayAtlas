from utils import *
from Kegg import *

def init_kegg_genome(recalc=False):
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
    init_kegg_genome()
    #  get all KEGG pathways and modules
    """
    kegg = KeggApi()
    pathways = kegg.get_all_pathways()
    modules = kegg.get_all_modules()
    """
    #  example get all variants in a pathway
    """
    pathway_obj = KeggNetwork('hsa01200', 'pathway')
    KeggNetwork.all_snvs(outpath='test.csv')
    """
