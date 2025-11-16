import os
import pickle
from typing import Optional, Dict

import metapredict as meta
import pandas as pd
from tqdm import tqdm

from definitions import KEGG_PATHWAY_MUTATIONS_PATH, DISORDERED_THRESHOLD
from utils import KeggGene


class DisorderPredict:
    def __init__(self, kegg):
        self.kegg = kegg

    def load_sequences_to_predict(self, sequences_pickle_file: str, recalc: bool=False) -> Optional[Dict[str, str]]:
        """
        Load pre-processed sequences from a pickle file or process them from scratch.

        Args:
            sequences_pickle_file (str): Path to pickle file.
            recalc (bool): Whether to force re-processing.

        Returns:
            Optional[Dict[str, str]]: Mapping from gene_id -> cleaned amino acid sequence.
        """
        if os.path.exists(sequences_pickle_file) and not recalc:
            print(f"Loading pre-processed sequences from cache: {sequences_pickle_file}")
            with open(sequences_pickle_file, 'rb') as f:
                return pickle.load(f)
        else:
            print("Processing sequences from scratch (or recalc=True)...")
            all_genes = list(self.kegg.get_all_genes().keys())
            sequences_to_predict = {}
            for gene_id in tqdm(all_genes, desc="Loading and cleaning gene sequences"):
                try:
                    gene = KeggGene(gene_id)
                    if gene.aa_seq and gene.coding_type == "CDS":
                        # Clean sequence based on metapredict V3 rules
                        processed_seq = gene.aa_seq.replace('B', 'N').replace('U', 'C').replace('X', 'G').replace('Z', 'Q')
                        processed_seq = processed_seq.replace(' ', '').replace('*', '').replace('-', '')
                        if processed_seq:
                            sequences_to_predict[gene_id] = processed_seq
                except Exception as e:
                    print(f"[Error] Loading gene {gene_id}: {e}")

            print(f"Saving {len(sequences_to_predict)} processed sequences to cache: {sequences_pickle_file}")
            with open(sequences_pickle_file, 'wb') as f:
                pickle.dump(sequences_to_predict, f)

            return sequences_to_predict

    @staticmethod
    def predict(sequences_to_predict: Dict[str, str]) -> None:
        """
        Predict disorder for the given sequences and merge results into SNV CSV files.

        Args:
            sequences_to_predict (Dict[str, str]): Mapping from gene_id -> amino acid sequence.
        """
        print(f"\n\n  Predicting disorder for {len(sequences_to_predict)} sequences...")
        try:
            disorder_predictions = meta.predict_disorder(sequences_to_predict)

            # Step 3: Loop through predictions, load corresponding SNV file, merge, and save.
            for gene_id, result_list in tqdm(disorder_predictions.items(),
                                             desc="Merging disorder scores into SNV files"):
                snv_file_path = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{gene_id}.csv")
                if not os.path.exists(snv_file_path):
                    continue
                snv_df = pd.read_csv(snv_file_path)

                scores = result_list[1]  # scores is a np array

                # Create a DataFrame from the disorder prediction list to act as a lookup table
                disorder_df = pd.DataFrame({
                    'amino_acid_index': range(len(scores)),
                    'disorder_score': scores
                })

                # Create the 'amino_acid_index' in the SNV DataFrame by dividing the nucleotide 'Start' position by 3.
                # This correctly maps the nucleotide position to the 0-indexed amino acid position.
                snv_df['amino_acid_index'] = snv_df['Start'] // 3

                # Merge the SNV data with the disorder scores
                merged_df = pd.merge(snv_df, disorder_df, on='amino_acid_index', how='left')

                merged_df['disorder_score'] = pd.to_numeric(merged_df['disorder_score'], errors='coerce')
                merged_df['disorder_score'] = merged_df['disorder_score'].fillna(0)  # NA -> ordered

                # Create the binary 'is_disordered' column
                merged_df['is_disordered'] = (merged_df['disorder_score'] > DISORDERED_THRESHOLD).astype(int)

                # Overwrite the original SNV file with the new, enriched data
                merged_df.to_csv(snv_file_path, index=False)

        except Exception as e:
            print(f"[Error] Failed during batch disorder prediction or file merging: {e}")
