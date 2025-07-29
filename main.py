from utils import *
from Kegg import *

from esm import pretrained



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


def genes_all_snvs(recalc=False):
    """Generate all SNVs for all KEGG genes and save them to files."""
    kegg = KeggApi()
    all_genes = list(kegg.get_all_genes().keys())
    if not recalc:
        all_genes = set(all_genes) - kegg_genes_in_dataset()
    all_genes = list(all_genes)
    print(f"Total genes to process: {len(all_genes)}")

    for gene_id in tqdm(all_genes, desc="Generating gene SNVs", unit="gene"):
        try:
            gene = KeggGene(gene_id)
            out_file = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{gene.kegg_id}.csv")
            if not os.path.exists(out_file):
                df = gene.all_snvs(outpath=out_file, index=True)
        except Exception as e:
            print(f"[Error] {gene_id}: {e}")


def pathways_and_modules_dict():
    """Create a dictionary for each KEGG pathway and module mapping gene IDs to their SNV files."""
    kegg = KeggApi()

    def process_collection(collection, network_type):
        for kegg_id in tqdm(collection, desc=f"Mapping {network_type}s", unit=network_type):
            try:
                gene_list = kegg.get_gene_list(kegg_id)
                gene_map = {}

                for gene_id in gene_list:
                    snv_file = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{gene_id}.csv")
                    if os.path.exists(snv_file):
                        gene_map[gene_id] = snv_file
                    else:
                        print(f"[Missing] SNV file not found for gene: {gene_id} in p/m: {kegg_id}")
                        # If SNV file is not found, we still want to keep the gene in the map
                        # The main reason of not having SNV file is that the len(gene)%3!=0
                        gene_map[gene_id] = None

                out_file = os.path.join(KEGG_PATHWAY_OBJECTS_PATH, f"{kegg_id}.pickle")
                pd.to_pickle(gene_map, out_file)
            except Exception as e:
                print(f"[Error] {kegg_id}: {e}")

    # Process pathways
    pathways = kegg.get_all_pathways()
    process_collection(pathways, 'pathway')

    # Process modules
    modules = kegg.get_all_modules()
    process_collection(modules, 'module')


def create_llr_scoring():
    """Create LLR scoring for all KEGG proteins using ESM1b model"""
    # Load ESM1b model
    model, alphabet = pretrained.load_model_and_alphabet(ESM1B_MODEL)
    calculator = WildtypeMarginalsCalculator(model=model, alphabet=alphabet)

    # Get KEGG genes
    kegg = KeggApi()
    all_genes = list(kegg.get_all_genes().keys())
    print(f"Total genes to process: {len(all_genes)}")

    # Get sequences
    seqs = kegg.gene_seq(all_genes, seq_type="aaseq")  # returns dict {gene_id: sequence}

    output_dir = "llr_scores"
    os.makedirs(output_dir, exist_ok=True)

    for gene_id in tqdm(seqs, desc="Generating LLR scoring", unit="gene"):
        seq = seqs[gene_id]
        if not seq or any(aa not in VALID_AA for aa in seq):
            print(f"Skipping {gene_id}: invalid or empty sequence")
            continue

        try:
            # Compute mutation scores
            scores = calculator.score_all_mutations(seq)
            output_path = os.path.join(output_dir, f"{gene_id.replace(':', '_')}_llr.csv")

            # Save scores to CSV
            with open(output_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Variant", "LLR_score"])
                for variant, score in scores.items():
                    writer.writerow([variant, score])

        except Exception as e:
            print(f"Error processing {gene_id}: {e}")






if __name__ == '__main__':
    # Lab Notebook:
    #    https://docs.google.com/document/d/1XR21LBpqW3q96BjExqsbH6JgEhV3Yc9sJuzG3abzUmY/edit?usp=sharing

    """
    # Step 1: Initialize KEGG genome (download all genes)
    print("Downloading/Loading KEGG genome...")
    init_kegg_genome()
    """

    """
    # step 2: create all_snvs for all genes [./data/kegg/pathways/snvs]
    print("Creating all SNVs in all genes...")
    genes_all_snvs()
    """

    # Step 3: build dict objects for each pathway from gene_id to snvs file [./data/kegg/pathways/objects]
    print("Mapping genes snvs to pathways and modules...")
    pathways_and_modules_dict()
    # Notice that some values in the dict may be None, which means that the SNV file is not available for that gene.

    # Step 4: create LLR scoring for all KEGG proteins using ESM1b model
    print("Creating LLR scoring for all KEGG proteins...")
    create_llr_scoring()
