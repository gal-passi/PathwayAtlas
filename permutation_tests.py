import pandas as pd
from statistical_models import *

class PermutatePathway:
    def __init__(self, pathways_files: List[str], cancer_scores_path: str, distance_metric: str,
                 n_permutations: int = BOOTSTRAP_SAMPLES, random_state: int = 42):

        for pathway_file in pathways_files:
            if not os.path.exists(pathway_file):
                print(f"Error: Pathway file does not exist: {pathway_file}")
                continue

        if not os.path.exists(cancer_scores_path):
            print(f"Error: Cancer results path does not exist: {cancer_scores_path}")
            return

        self.cancer_results_path = cancer_scores_path
        self.pathways_files = pathways_files
        self.distance_metric = distance_metric
        self.n_permutations = n_permutations
        self.random_state = random_state
        self.gmm_fitter = GMM()

    def run_test(self):
        all_cancer_score_files = glob.glob(pjoin(self.cancer_results_path, "*.csv"))
        if not all_cancer_score_files:
            print(f"No result CSV files found in the specified directory: {self.cancer_results_path}")
            return

        cancers_df_dict = {}

        for cancer_scores_file in all_cancer_score_files:
            df = pd.read_csv(cancer_scores_file)
            cancers_df_dict[cancer_scores_file] = df

            for col in ['p_value', 'q_value'] + [f'significant_{alpha}' for alpha in P_VALUE_THRESHOLDS]:
                if col not in df.columns:
                    df[col] = np.nan

        for pathway_file in self.pathways_files:
            pathway_name = os.path.splitext(os.path.basename(pathway_file))[0]
            pathway_df = pd.read_csv(pathway_file)
            p_values = []

            print(f"======= Running permutation test for pathway: {pathway_name} =======")

            for cancer_scores_file, cancer_df in cancers_df_dict.items():

                print(f"  Processing cancer file: {cancer_scores_file}")

                # Filter the cancer DataFrame for the current pathway
                pathway_cancer_df = cancer_df[cancer_df['pathway_name'] == pathway_name]

                if pathway_cancer_df.empty:
                    print(f"No data for pathway {pathway_name} in cancer file {cancer_scores_file}. Skipping.")
                    continue

                num_samples_dict = self._get_num_samples_per_protein(pathway_df, cancer_df)
                distances = self._bootstrap_pathway(pathway_cancer_df, num_samples_dict)

                observed_distance = pathway_cancer_df['distance'].iloc[0]

                p_value = self._calculate_p_value(observed_distance, distances)
                cancer_df.loc[pathway_cancer_df.index, 'p_value'] = p_value
                p_values.append(p_value)

                print(f"    Observed distance: {observed_distance}, p-value: {p_value}")

            # FDR correction
            reject_list, q_values = self._calculate_q_value(p_values)

            for i, cancer_scores_file in enumerate(cancers_df_dict.keys()):
                cancer_df = cancers_df_dict[cancer_scores_file]
                pathway_indices = cancer_df[cancer_df['pathway_name'] == pathway_name].index
                if not pathway_indices.empty:
                    cancer_df.loc[pathway_indices, 'q_value'] = q_values[i]
                    for alpha in P_VALUE_THRESHOLDS:
                        significant_col = f'significant_{alpha}'
                        cancer_df.loc[pathway_indices, significant_col] = reject_list[i]

    def _bootstrap_pathway(self, bg_scores_df: pd.DataFrame, num_samples_dict: Dict[str, int]) -> List[float]:
        """
        Bootstrapping logic for a specific pathway.
        """
        distances = []

        # --- 1. PREPARE REFERENCE DISTRIBUTION (ONCE) ---
        ref_scores_dict = self.get_pathway_scores_background(bg_scores_df)

        # Create weighted histogram (PSSM)
        ref_hist, bin_edges = CancerPathwayScoring.create_joint_distribution(
            ref_scores_dict, MICHAL_HN1_PSSM, use_pssm=True
        )

        # If KL, we need to fit the Reference GMM once
        ref_gmm = None
        if self.distance_metric == "kl_divergence":
            gmm_fitter = GMM()
            ref_gmm, _ = gmm_fitter.GMM_the_distribution(ref_hist, bin_edges)
            if ref_gmm is None:
                print("    Error: Could not fit GMM to background reference.")
                return []

        # --- 2. BOOTSTRAP LOOP ---
        for i in range(self.n_permutations):
            # A. Sample from the dataframe

            sampled_df = self._sample_pathway_by_proteins(bg_scores_df, num_samples_dict)

            # B. Create Sampled Distribution (Weighted by PSSM same as observed data logic)
            sampled_scores_dict = self.get_pathway_scores_background(sampled_df)
            sampled_hist, bin_edges = CancerPathwayScoring.create_joint_distribution(
                sampled_scores_dict, MICHAL_HN1_PSSM, use_pssm=True
            )

            # Check for empty distributions (prevent crash)
            if np.sum(sampled_hist) == 0:
                continue

            # C. Calculate Metric
            dist = 0.0

            if self.distance_metric == "kl_divergence":
                # Fit GMM for sample
                sampled_gmm, _ = self.gmm_fitter.GMM_the_distribution(sampled_hist, bin_edges)
                if sampled_gmm is not None:
                    dist = DistributionDistances.kl_divergence_from_gmms(
                        ref_gmm, sampled_gmm, n_samples_mc=RANDOM_SAMPLE_NUM
                    )
                else:
                    continue  # Skip if GMM fails

            elif self.distance_metric == "wasserstein":
                # Direct Histogram calculation
                dist = DistributionDistances.wasserstein_from_hist(ref_hist, sampled_hist, bin_edges)

            elif self.distance_metric in ["dw_distance", "directional_wasserstein"]:
                # Returns (distance, shift). We bootstrap the MAGNITUDE (distance)
                # to see if the deviation is significant.
                _, dw_shift = DistributionDistances.directional_wasserstein_from_hist(
                    ref_hist, sampled_hist, bin_edges
                )
                dist = dw_shift

            distances.append(dist)

            if i % 100 == 0 and i > 0:
                print(f"    Permutation {i}/{self.n_permutations}...")

        return distances

    @staticmethod
    def _get_num_samples_per_protein(pathway_df: pd.DataFrame, cancer_df: pd.DataFrame) -> Dict[str, int]:
        """
        Counts the number of unique mutations per protein in the cancer DataFrame
        that are present in the pathway DataFrame.
        Returns
        -------
        dict
            A dictionary mapping protein names to the count of unique mutations.
        """
        if pathway_df.empty or cancer_df.empty:
            print("Warning: Empty DataFrame provided.")
            return {}

        samples_dict = {}
        protein_names = pathway_df[PROTEIN_NAME_COL].unique()

        for protein in protein_names:
            protein_mutations = cancer_df[cancer_df[PROTEIN_NAME_COL] == protein]
            num_unique_mutations = protein_mutations[VARIANT_COL].nunique()
            samples_dict[protein] = num_unique_mutations

    def _sample_pathway_by_proteins(self, bg_scores_df: pd.DataFrame, num_samples_dict: Dict[str, int]) -> pd.DataFrame:
        """
        Samples the background scores DataFrame based on the number of samples
        required per protein.
        Returns
        -------
        pd.DataFrame
            A DataFrame containing the sampled rows.
        """
        sampled_rows = []

        for protein, num_samples in num_samples_dict.items():
            protein_bg_scores = bg_scores_df[bg_scores_df[PROTEIN_NAME_COL] == protein]
            if protein_bg_scores.empty:
                continue

            sampled_protein_scores = protein_bg_scores.sample(
                n=num_samples, replace=True, random_state=self.random_state
            )
            sampled_rows.append(sampled_protein_scores)

        if sampled_rows:
            return pd.concat(sampled_rows, ignore_index=True)
        else:
            return pd.DataFrame(columns=bg_scores_df.columns)

    @staticmethod
    def _calculate_p_value(observed_distance, bootstrap_distances):
        """
        Calculates the p-value based on the observed distance and the distribution
        of distances from permutations.

        Args:
            observed_distance (float): The distance calculated from the actual data.
            bootstrap_distances (list): List of distances from permuted datasets.

        Returns:
            float: The p-value indicating the significance of the observed distance.
        """
        if not bootstrap_distances:
            return 1.0

        # Convert to numpy array for efficiency
        dist_array = np.array(bootstrap_distances)

        # P-value = (Number of random dists >= observed dist) / Total permutations
        # We add 1 to numerator and denominator for pseudo-count smoothing (prevents p=0)
        n_extreme = np.sum(dist_array >= observed_distance)
        return (n_extreme + 1) / (len(dist_array) + 1)

    @staticmethod
    def _calculate_q_value(p_values: list, alpha: float = 0.05) -> tuple:

        """
        Applies Benjamini-Hochberg FDR correction to a list of p-values.

        Args:
            p_values (list): List of p-values to correct.
            alpha (float): Significance level for FDR correction.

        Returns:
            list: List of FDR-corrected p-values.
        """
        corrected_p_values = []
        reject = []

        for p in p_values:
            if p < 0 or p > 1 or np.isnan(p):
                corrected_p_values.append(np.nan)
                reject.append(False)
            else:
                q = false_discovery_control(p, method='bh')
                corrected_p_values.append(q)
                r = q <= alpha
                reject.append(r)

        return reject, corrected_p_values

    @staticmethod
    def get_pathway_scores_background(pathway_df) -> dict:
        """
        Collects mutation scores from the pathway's genes, categorized by
        mutation type (e.g., 'A>C') and then by score type.

        Returns
        -------
        dict
            A nested dictionary: {mut_type: {score_type: [scores]}}.
            Example: {"A>C": {"esm_log_probs": [0.1, 0.2], ...}, ...}
        """
        if pathway_df.empty:
            return {}

        background_scores = defaultdict(lambda: defaultdict(list))

        # Check for essential columns in the pre-aggregated CSV
        required_cols = {"Ref", "Alt"}
        if not required_cols.issubset(pathway_df.columns):
            print("Warning: 'Ref' or 'Alt' columns missing in pathway scores file")
            return {}

        score_types = ['esm_log_probs', 'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob']

        # Efficiently iterate over the single, preloaded DataFrame
        for _, row in pathway_df.iterrows():
            mut_type = f"{str(row['Ref']).upper()}>{str(row['Alt']).upper()}"

            for score_type in score_types:
                # Check if score exists and is not NaN
                if score_type in row and pd.notna(row[score_type]):
                    background_scores[mut_type][score_type].append(float(row[score_type]))

        # Convert default dicts to regular dicts for a clean return value
        return {mut: dict(scores) for mut, scores in background_scores.items()}
