from typing import List, Optional, Dict

from definitions import *
import pandas as pd
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
import pickle
from scipy.stats import wasserstein_distance


class GMM:
    def GMM_the_distribution(self, joint_distribution, bin_edges,
                             n_components=GMM_COMPONENTS, n_samples=RANDOM_SAMPLE_NUM, random_state=None):
        """
        Fit a Gaussian Mixture Model (GMM) to a distribution defined by a histogram.
        Returns the fitted GMM and its BIC.
        """
        # 1. Reconstruct the dataset from the histogram
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        counts = (joint_distribution * n_samples).astype(int)

        if np.sum(counts) == 0:
            return None

        reconstructed_data = np.repeat(bin_centers, counts).reshape(-1, 1)

        # 2. fit model
        gmm = GaussianMixture(
            n_components=n_components,
            random_state=random_state,
            covariance_type='full'
        )
        gmm.fit(reconstructed_data)

        bic = gmm.bic(reconstructed_data)
        return gmm, bic

    def GMM_bic_curve(self, joint_distribution, bin_edges, max_components=10,
                      n_samples=RANDOM_SAMPLE_NUM, random_state=None, filename='bic_curve'):
        """
        Plot the BIC curve for different numbers of Gaussian components and save it.

        Parameters
        ----------
        joint_distribution : np.ndarray
            Histogram values of the distribution.
        bin_edges : np.ndarray
            Histogram bin edges.
        max_components : int
            Maximum number of GMM components to test.
        n_samples : int
            Number of samples to reconstruct the dataset from the histogram.
        random_state : int, optional
            Random state for reproducibility.
        filename : str
            The base name for the saved plot file (e.g., 'my_bic_plot').
            The extension '.png' will be added automatically.
        """
        # 1. Reconstruct dataset
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        counts = (joint_distribution * n_samples).astype(int)

        if np.sum(counts) == 0:
            print("Empty distribution, skipping GMM fitting.")
            return

        reconstructed_data = np.repeat(bin_centers, counts).reshape(-1, 1)

        # 2. Fit GMMs with different n_components and record BIC
        bics = []
        n_components_range = range(1, max_components + 1)

        for n_components in n_components_range:
            gmm = GaussianMixture(
                n_components=n_components,
                random_state=random_state,
                covariance_type='full'
            )
            gmm.fit(reconstructed_data)
            bics.append(gmm.bic(reconstructed_data))

        # 3. Plot the BIC curve
        plt.figure(figsize=(6, 4))
        plt.plot(n_components_range, bics, marker='o')
        plt.xlabel("Number of components")
        plt.ylabel("BIC score")
        plt.title("GMM BIC Curve")
        plt.grid(True)

        # Construct the full save path
        save_dir = "bic_for_selecting_num_of_gmm"
        full_save_path = os.path.join(save_dir, f"{filename}.png")

        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            print(f"Created directory: {save_dir}")

        # Save the figure
        try:
            plt.savefig(full_save_path, bbox_inches='tight', dpi=300)
            print(f"BIC curve saved to: {full_save_path}")
        except Exception as e:
            print(f"Error saving BIC curve to {full_save_path}: {e}")

    def plot_1D_GMM(self, joint_distribution, bin_edges, gmm, ax=None):
        """
        Plots the histogram of the joint distribution and overlays the
        probability density function of the fitted Gaussian Mixture Model (GMM).

        This visualization helps to assess the quality of the GMM fit to the
        underlying distribution. It also plots the individual Gaussian components
        that make up the mixture model.

        Parameters
        ----------
        joint_distribution : np.ndarray
            The values (heights) of the histogram representing the distribution.
        bin_edges : np.ndarray
            The edges of the histogram bins.
        gmm : sklearn.mixture.GaussianMixture
            The fitted Gaussian Mixture Model object to be plotted.
        ax : matplotlib.axes.Axes, optional
            The matplotlib axes object on which to plot. If None, a new figure
            and axes will be created, by default None.

        Returns
        -------
        matplotlib.axes.Axes
            The axes object containing the plot.
        """
        # If no axes are provided, create a new figure and axes
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        # 1. Plot the original distribution histogram
        # Calculate bin widths and centers for plotting the histogram
        bin_widths = np.diff(bin_edges)
        bin_centers = bin_edges[:-1] + bin_widths / 2

        # Plot the histogram using a bar plot. The height is the probability density.
        # The area of the histogram should sum to 1.
        ax.bar(bin_centers, joint_distribution, width=bin_widths, color='gray', alpha=0.6,
               label='Original Distribution')

        # 2. Plot the overall GMM probability density function
        # Create a fine grid of x-values for a smooth plot of the PDF
        x_grid = np.linspace(bin_edges[0], bin_edges[-1], 1000).reshape(-1, 1)

        # The score_samples method returns the log-likelihood of each sample
        log_pdf = gmm.score_samples(x_grid)
        pdf = np.exp(log_pdf)

        # Plot the total GMM PDF
        ax.plot(x_grid, pdf, color='red', linewidth=2, label='GMM Fit')

        # 3. Plot the individual Gaussian components of the GMM
        # This shows how the GMM is composed
        for i in range(gmm.n_components):
            # Extract the parameters for the i-th component
            mean = gmm.means_[i][0]
            # In scikit-learn's GMM, covariance is stored as a matrix
            std_dev = np.sqrt(gmm.covariances_[i][0][0])
            weight = gmm.weights_[i]

            # Calculate the PDF for this individual component, scaled by its weight
            component_pdf = weight * norm.pdf(x_grid.flatten(), mean, std_dev)

            # Plot the component PDF with a dashed line
            ax.plot(x_grid, component_pdf, linestyle='--', label=f'Component {i + 1}')

        # 4. Final plot formatting
        ax.set_title('Gaussian Mixture Model Fit to Distribution')
        ax.set_xlabel('Value')
        ax.set_ylabel('Probability Density')
        ax.legend()
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.show()

    # TODO use this in the results functions
    @staticmethod
    def save_distribution(pathway_id, joint_distribution, bin_edges,
                          gmm_model=None, output_dir=DISTRIBUTIONS_PATH):
        """
        Saves the joint distribution, bin edges, and the fitted GMM model to a file using pickle.

        Parameters
        ----------
        pathway_id : str
            An identifier for the pathway (e.g., "hsa00010"), used for the filename.
        joint_distribution : np.ndarray
            The values (heights) of the histogram representing the distribution.
        bin_edges : np.ndarray
            The edges of the histogram bins.
        gmm_model : sklearn.mixture.GaussianMixture, optional
            The fitted Gaussian Mixture Model object. If None, it will not be saved.
        output_dir : str, optional

        Returns
        -------
        None
        """
        # Ensure the output directory exists. If not, create it.
        os.makedirs(output_dir, exist_ok=True)

        # Prepare the data to be saved in a dictionary for clarity and structure.
        data_to_save = {
            'pathway_id': pathway_id,
            'joint_distribution': joint_distribution,
            'bin_edges': bin_edges,
            'gmm_model': gmm_model
        }

        file_path = os.path.join(output_dir, f"{pathway_id}_distribution.pickle")

        try:
            with open(file_path, 'wb') as f:
                pickle.dump(data_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Successfully saved distribution for {pathway_id} to {file_path}")
        except Exception as e:
            print(f"Error saving distribution for {pathway_id}: {e}")




class DistributionDistances:
    """
    A utility class to calculate distances between two data distributions.

    Contains methods for:
    - Chi-Squared distance from histograms.
    - Kullback-Leibler (KL) Divergence from fitted Gaussian Mixture Models (GMMs).
    - 1-Wasserstein distance (Earth Mover's Distance)
    """

    @staticmethod
    def chi_squared_from_hist(counts1: np.ndarray, counts2: np.ndarray) -> float:
        """
        Calculates the Chi-Squared distance between two pre-computed histograms.
        This method assumes that both histograms were created using the exact same
        set of bins.
        Args:
            counts1 (np.ndarray): An array of bin counts for the first histogram.
            counts2 (np.ndarray): An array of bin counts for the second histogram.
        Returns:
            float: The calculated Chi-Squared distance. A smaller value indicates
                   more similar distributions.
        Raises:
            ValueError: If the input count arrays do not have the same length.
        """
        if len(counts1) != len(counts2):
            raise ValueError("Input histogram count arrays must have the same length.")

        # Ensure counts are float to avoid integer division issues
        counts1 = counts1.astype(float)
        counts2 = counts2.astype(float)

        numerator = (counts1 - counts2) ** 2
        denominator = counts1 + counts2

        # Handle bins where the denominator is zero (i.e., the bin is empty for both)
        # The contribution of these bins to the distance is 0.
        chi_squared_terms = np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),  # Where denominator is 0, the output is 0
            where=(denominator != 0)  # Condition for the division
        )

        return np.sum(chi_squared_terms)

    @staticmethod
    def kl_divergence_from_gmms(gmm_p: GaussianMixture, gmm_q: GaussianMixture,
                                n_samples_mc: int = RANDOM_SAMPLE_NUM) -> float:
        """
        Calculates KL-Divergence KL(P || Q) between two pre-fitted GMMs.

        This is the recommended method when you already have fitted GMMs.

        Args:
            gmm_p (GaussianMixture): The fitted GMM for distribution P.
            gmm_q (GaussianMixture): The fitted GMM for distribution Q.
            n_samples_mc (int): Number of samples for Monte Carlo approximation.

        Returns:
            float: The estimated KL-Divergence KL(P || Q).
        """
        # 1. Draw samples from the first GMM (P)
        samples, _ = gmm_p.sample(n_samples_mc)

        # 2. Calculate log-likelihood of samples under both models
        log_p_p = gmm_p.score_samples(samples)
        log_q_p = gmm_q.score_samples(samples)

        # 3. The KL divergence is the expectation E[log(P(x)) - log(Q(x))]
        return np.mean(log_p_p - log_q_p)

    @staticmethod
    def wasserstein_from_hist(counts1: np.ndarray, counts2: np.ndarray, bin_edges: np.ndarray) -> float:
        """
        Calculates the 1-Wasserstein distance (Earth Mover's Distance) between
        two pre-computed histograms.

        This metric represents the minimum "cost" to transform one distribution
        into the other, where cost is the amount of probability mass moved
        multiplied by the distance it is moved.

        Args:
            counts1 (np.ndarray): An array of bin densities for the first histogram.
            counts2 (np.ndarray): An array of bin densities for the second histogram.
            bin_edges (np.ndarray): The edges of the histogram bins.

        Returns:
            float: The calculated 1-Wasserstein distance.
        """
        if len(counts1) != len(counts2):
            raise ValueError("Input histogram count arrays must have the same length.")

        # Calculate the center of each bin to represent the location of the mass
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # SciPy's function calculates the distance between two 1D distributions.
        # It takes the locations (values) and their respective weights (probabilities).
        return wasserstein_distance(u_values=bin_centers, v_values=bin_centers,
                                    u_weights=counts1, v_weights=counts2)





class CancerPathwayScoring:
    """
    Analyzes and compares mutation scores between a KEGG pathway (background model)
    and a specific cancer type's mutation dataset.
    """

    def __init__(self, pathway_dict_path: str, pathway_scores_csv_path: str, cancer_csv_path: str):
        """
        Initializes the scoring object by loading pathway and cancer data.

        Parameters
        ----------
        pathway_dict_path : str
            File path to the pickled pathway dictionary. Used to get the list of gene IDs.
        pathway_scores_csv_path : str
            File path to the pre-aggregated CSV with all background scores for the pathway.
        cancer_csv_path : str
            File path to the CSV with scored cancer mutations.
        """
        self.pathway_dict_path = pathway_dict_path
        self.pathway_scores_csv_path = pathway_scores_csv_path
        self.cancer_csv_path = cancer_csv_path
        self.pssm = MICHAL_HN1_PSSM  # Using the predefined PSSM for weighting
        self.gmm_fitter = GMM()  # GMM helper instance for fitting models

        # --- MODIFICATION 1: Load both pathway definition and pathway scores ---
        # Load the pickle file to get the definitive list of gene IDs
        try:
            with open(self.pathway_dict_path, 'rb') as f:
                self.pathway_dict = pickle.load(f)
        except (FileNotFoundError, EOFError) as e:
            print(f"[Error] Could not load pathway dictionary {self.pathway_dict_path}: {e}")
            self.pathway_dict = {}

        # Load the pre-aggregated CSV with all background scores
        try:
            self.pathway_df = pd.read_csv(self.pathway_scores_csv_path)
        except FileNotFoundError:
            print(f"[Error] Pathway scores CSV not found: {self.pathway_scores_csv_path}")
            self.pathway_df = pd.DataFrame()

        # Load cancer data (no change here)
        try:
            self.cancer_df = pd.read_csv(self.cancer_csv_path)
        except FileNotFoundError:
            print(f"[Error] Cancer CSV file not found: {self.cancer_csv_path}")
            self.cancer_df = pd.DataFrame()

    def get_pathway_genes_id(self) -> set:
        """Returns a set of KEGG gene IDs in the pathway."""
        return set(self.pathway_dict.keys())

    def get_pathway_scores_background(self) -> dict:
        """
        Collects mutation scores from the pathway's genes, categorized by
        mutation type (e.g., 'A>C') and then by score type.

        Returns
        -------
        dict
            A nested dictionary: {mut_type: {score_type: [scores]}}.
            Example: {"A>C": {"esm_log_probs": [0.1, 0.2], ...}, ...}
        """
        if self.pathway_df.empty:
            return {}

        background_scores = defaultdict(lambda: defaultdict(list))

        # Check for essential columns in the pre-aggregated CSV
        required_cols = {"Ref", "Alt"}
        if not required_cols.issubset(self.pathway_df.columns):
            print(f"Warning: 'Ref' or 'Alt' columns missing in pathway scores file: {self.pathway_scores_csv_path}")
            return {}

        score_types = ['esm_log_probs', 'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob']

        # Efficiently iterate over the single, pre-loaded DataFrame
        for _, row in self.pathway_df.iterrows():
            mut_type = f"{str(row['Ref']).upper()}>{str(row['Alt']).upper()}"

            for score_type in score_types:
                # Check if score exists and is not NaN
                if score_type in row and pd.notna(row[score_type]):
                    background_scores[mut_type][score_type].append(float(row[score_type]))

        # Convert defaultdicts to regular dicts for a clean return value
        return {mut: dict(scores) for mut, scores in background_scores.items()}

    def get_cancer_pathway_scores(self) -> tuple[dict, int]:
        """
        Filters the cancer dataset for mutations in the pathway's genes and
        collects their scores.

        Returns
        -------
        tuple[dict, int]
            A tuple containing:
            - A nested dictionary of scores: {mut_type: {score_type: [scores]}}.
            - An integer count of the relevant mutations found.
        """
        cancer_scores = defaultdict(lambda: defaultdict(list))
        pathway_gene_ids = self.get_pathway_genes_id()

        if self.cancer_df.empty or not pathway_gene_ids:
            return {}, 0  # Return count of 0

        required_cols = {'KeggId', 'Ref', 'Alt'}
        if not required_cols.issubset(self.cancer_df.columns):
            return {}, 0  # Return count of 0

        relevant_df = self.cancer_df.dropna(subset=['KeggId', 'Ref', 'Alt']).copy()

        def is_in_pathway(kegg_id_cell: str) -> bool:
            return any(gid in pathway_gene_ids for gid in str(kegg_id_cell).split(','))

        mask = relevant_df['KeggId'].apply(is_in_pathway)
        pathway_cancer_mutations = relevant_df[mask]

        mutation_count = len(pathway_cancer_mutations)

        if pathway_cancer_mutations.empty:
            return {}, 0  # Return count of 0

        for _, row in pathway_cancer_mutations.iterrows():
            mut_type = f"{str(row['Ref']).upper()}>{str(row['Alt']).upper()}"
            for score_type in ['esm_log_probs', 'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob']:
                if score_type in row and pd.notna(row[score_type]):
                    cancer_scores[mut_type][score_type].append(float(row[score_type]))

        return {mut: dict(scores) for mut, scores in cancer_scores.items()}, mutation_count

    @staticmethod
    def create_joint_distribution(scores_dict, pssm,
                                  score_type="clinvar_reg_dis_ordered_prob",
                                  bins=BINS_NUM, range=(0, 1), use_pssm: bool = True):
        """
        Creates a joint distribution from a scores dictionary.
        If use_pssm is True, it creates a PSSM-weighted distribution.
        Otherwise, it creates a simple, unweighted distribution.
        """
        if isinstance(bins, int):
            bin_edges = np.linspace(range[0], range[1], bins + 1)
            num_bins = bins
        else:  # bins is already an array of edges
            bin_edges = bins
            num_bins = len(bin_edges) - 1

        if use_pssm:
            # PSSM-weighted logic for the background distribution
            joint_distribution = np.zeros(num_bins)
            total_weight_used = 0.0

            for mut_type, weight in pssm.items():
                if mut_type in scores_dict and score_type in scores_dict[mut_type]:
                    valid_scores = [s for s in scores_dict[mut_type][score_type] if s is not None and np.isfinite(s)]
                    if not valid_scores:
                        continue
                    hist, _ = np.histogram(valid_scores, bins=bin_edges, density=True)
                    joint_distribution += hist * weight
                    total_weight_used += weight
            if total_weight_used > 0:
                joint_distribution /= total_weight_used
        else:
            # Unweighted logic for the observed cancer distribution
            all_scores = []
            for mut_type in scores_dict:
                if score_type in scores_dict[mut_type]:
                    valid_scores = [s for s in scores_dict[mut_type][score_type] if s is not None and np.isfinite(s)]
                    all_scores.extend(valid_scores)

            if not all_scores:
                return np.zeros(num_bins), bin_edges

            joint_distribution, _ = np.histogram(all_scores, bins=bin_edges, density=True)

        return joint_distribution, bin_edges

    def get_bins_of_distributions_ready(self, background_scores: dict, cancer_scores: dict,
                                        score_type: str = "clinvar_reg_dis_ordered_prob"):
        """
        A helper method to generate the background and cancer distributions.
        This centralizes the distribution creation logic to be used by multiple
        distance calculation methods.

        Returns
        -------
        tuple
            (background_distribution, bin_edges, cancer_distribution)
        """
        # 1. Create a PSSM-weighted joint distribution for the background model
        bg_dist, bin_edges = self.create_joint_distribution(
            background_scores, self.pssm, score_type=score_type, use_pssm=True
        )

        # 2. Create a simple, unweighted distribution for the cancer scores
        cancer_dist, _ = self.create_joint_distribution(
            cancer_scores, self.pssm, score_type=score_type, bins=bin_edges, use_pssm=False
        )

        return bg_dist, bin_edges, cancer_dist

    def calculate_distance_gmm_kl_d(self, background_scores: dict, cancer_scores: dict,
                                    score_type: str = "clinvar_reg_dis_ordered_prob") -> float:
        """
        Calculates KL-Divergence between background and cancer score distributions
        by modeling them with GMMs.
        """
        # 1. Get the distributions from the helper method
        bg_dist, bin_edges, cancer_dist = self.get_bins_of_distributions_ready(
            background_scores, cancer_scores, score_type
        )

        if np.sum(bg_dist) == 0 or np.sum(cancer_dist) == 0:
            print(f"Warning: Empty distribution for '{score_type}'. Cannot calculate distance.")
            return None

        # 2. Fit GMM to each distribution
        gmm_bg, _ = self.gmm_fitter.GMM_the_distribution(bg_dist, bin_edges)
        gmm_cancer, _ = self.gmm_fitter.GMM_the_distribution(cancer_dist, bin_edges)

        if gmm_bg is None or gmm_cancer is None:
            print(f"Warning: GMM fitting failed for '{score_type}'.")
            return None

        # 3. Calculate KL-Divergence between the fitted models
        kl_divergence = DistributionDistances.kl_divergence_from_gmms(gmm_bg, gmm_cancer)

        return kl_divergence

    def calculate_distance_wasserstein(self, background_scores: dict, cancer_scores: dict,
                                       score_type: str = "clinvar_reg_dis_ordered_prob") -> float:
        """
        Calculates the 1-Wasserstein distance between the background and cancer
        score distributions directly from their histograms.
        """
        # 1. Get the distributions from the helper method
        bg_dist, bin_edges, cancer_dist = self.get_bins_of_distributions_ready(
            background_scores, cancer_scores, score_type
        )

        if np.sum(bg_dist) == 0 or np.sum(cancer_dist) == 0:
            print(f"Warning: Empty distribution for '{score_type}'. Cannot calculate distance.")
            return None

        # 2. Calculate Wasserstein distance directly from the histogram densities
        w_distance = DistributionDistances.wasserstein_from_hist(bg_dist, cancer_dist, bin_edges)

        return w_distance




class PathwayAtlasResults:
    """
    This class orchestrates the calculation of distance scores for all pairs of pathway-cancer data.
    """

    def __init__(self, cancer_csvs_path: str = CANCER_CSVS_MUTATIONS,
                 pathways_dicts_path: str = KEGG_PATHWAY_OBJECTS_PATH,
                 pathways_scores_path: str = PATHWAY_SCORES_PATH):
        """
        Initializes the PathwayAtlasResults with paths to the data directories.

        Args:
            cancer_csvs_path (str): Path to the directory containing cancer mutation CSV files.
            pathways_dicts_path (str): Path to the directory containing pathway data files (e.g., pickles).
        """
        if not os.path.isdir(cancer_csvs_path):
            raise FileNotFoundError(f"The specified cancer CSVs path does not exist: {cancer_csvs_path}")
        if not os.path.isdir(pathways_dicts_path):
            raise FileNotFoundError(f"The specified pathways dictionaries path does not exist: {pathways_dicts_path}")

        self.__cancer_csvs_path = cancer_csvs_path
        self.__pathways_dicts_path = pathways_dicts_path
        self.__pathways_scores_path = pathways_scores_path
        self.results = {}
        print("PathwayAtlasResults initialized successfully.")
        print(f"Cancer data will be read from: {self.__cancer_csvs_path}")
        print(f"Pathway data will be read from: {self.__pathways_dicts_path}\n\n")

    def _calculate_distance(self, analyzer, scoring_system: str, background_scores: Dict,
                            cancer_scores: Dict, score_to_analyze: str) -> Optional[float]:
        """
        Calculates the distance based on the specified scoring system.

        Args:
            analyzer: An instance of CancerPathwayScoring.
            scoring_system (str): The scoring system to use ('kl_divergence' or 'wasserstein').
            background_scores (Dict): A dictionary of background scores for genes in the pathway.
            cancer_scores (Dict): A dictionary of cancer-specific scores for genes in the pathway.
            score_to_analyze (str): The key for the score to extract from the score dictionaries.

        Returns:
            Optional[float]: The calculated distance, or None if the scoring system is invalid.
        """
        distance_calculators = {
            "kl_divergence": analyzer.calculate_distance_gmm_kl_d,
            "wasserstein": analyzer.calculate_distance_wasserstein,
        }
        calculator = distance_calculators.get(scoring_system)
        if calculator:
            return calculator(background_scores, cancer_scores, score_type=score_to_analyze)
        else:
            print(f"Error: Invalid scoring system specified: '{scoring_system}'. Cannot calculate distance.")
            return None

    def calculate_distances_single_cancer(self, cancer_file_path: str, score_to_analyze: str,
                                          scoring_system: str) -> Dict[str, float]:
        """
        For each pathway, calculates the distance to a single cancer dataset.

        This function iterates through all available pathway definition files (pickles),
        finds their corresponding aggregated score files (CSVs), and computes a
        distance score against the provided cancer mutation data.

        Args:
            cancer_file_path (str): The file path for a single cancer's mutation data.
            score_to_analyze (str): The specific score to be analyzed (e.g., 'clinvar_reg_dis_ordered_prob').
            scoring_system (str): The distance metric to use, e.g., 'kl_divergence' or 'wasserstein'.

        Returns:
            Dict[str, float]: A dictionary where keys are pathway file names and values are the calculated distances.
        """
        if not os.path.isfile(cancer_file_path):
            print(f"Error: Cancer file not found at path: {cancer_file_path}")
            return {}

        cancer_name = os.path.basename(cancer_file_path)
        print(f"\n--- Starting Analysis for Cancer: {cancer_name} ---")
        self.results[cancer_name] = {}

        # Get the list of all pathway definition files (pickles)
        pathway_dict_files = sorted(os.listdir(self.__pathways_dicts_path))

        for i, pathway_dict_filename in enumerate(pathway_dict_files):
            # Construct the full path to the pathway's pickle file
            pathway_dict_filepath = os.path.join(self.__pathways_dicts_path, pathway_dict_filename)

            # From the pickle filename, derive the corresponding scores CSV filename
            base_name = os.path.splitext(pathway_dict_filename)[0]
            pathway_scores_filename = f"{base_name}.csv"
            pathway_scores_filepath = os.path.join(self.__pathways_scores_path, pathway_scores_filename)

            print(f"\n({i + 1}/{len(pathway_dict_files)}) Analyzing Pathway: {base_name}")

            # Safety check: ensure the required scores CSV file actually exists
            if not os.path.isfile(pathway_scores_filepath):
                print(f"--> Warning: Scores CSV '{pathway_scores_filename}' not found [empty pathway probably]. Skipping.")
                continue

            try:
                # Instantiate the analyzer with the three required file paths
                analyzer = CancerPathwayScoring(
                    pathway_dict_path=pathway_dict_filepath,
                    pathway_scores_csv_path=pathway_scores_filepath,
                    cancer_csv_path=cancer_file_path
                )

                # Get the score dictionaries. The background one is now loaded very fast.
                background_scores = analyzer.get_pathway_scores_background()
                cancer_scores, mutation_count = analyzer.get_cancer_pathway_scores()

                # Proceed only if both datasets contain relevant data
                if background_scores and cancer_scores:
                    distance = self._calculate_distance(analyzer, scoring_system,
                                                        background_scores, cancer_scores,
                                                        score_to_analyze)
                    if distance is not None:
                        # Store the result: { "hsa04010.pickle": ... }
                        # Store a dictionary with both distance and count
                        self.results[cancer_name][pathway_dict_filename] = {
                            'distance': distance,
                            'n': mutation_count
                        }
                        print(f"--> Success! Distance ({scoring_system}): {distance:.4f}   (n={mutation_count})")
                else:
                    # Provide a clear reason for skipping
                    print(f"--> Analysis skipped.")
                    if not background_scores:
                        print("    Reason: No background scores could be loaded from the pathway files.")
                    if not cancer_scores:
                        print("    Reason: No mutations from this cancer were found in the pathway's genes.")

            except Exception as e:
                print(f"--> An unexpected error occurred during analysis of {cancer_name}-{base_name}: {e}")

        print(f"\n--- Analysis Complete for Cancer: {cancer_name} ---")
        return self.results[cancer_name]

    # this might take forever to run, split to sarray jobs
    # def run_full_analysis(self, score_to_analyze: str, scoring_system: str,
    #                       results_path: str = RESULTS_PATH) -> Dict[str, Dict[str, float]]:
    #     """
    #     Runs the entire analysis by iterating through all cancer files,
    #     calculating distances against all pathways, and saving backup results incrementally.
    #
    #     Args:
    #         score_to_analyze (str): The score type to use for the analysis
    #                                 (e.g., 'clinvar_reg_dis_ordered_prob').
    #         scoring_system (str): The distance metric to use
    #                               (e.g., 'kl_divergence', 'wasserstein').
    #         results_path (str): The directory where backup and final results will be saved.
    #
    #     Returns:
    #         Dict[str, Dict[str, float]]: A nested dictionary containing the complete results.
    #         The structure is: {cancer_name: {pathway_file: distance_score}}.
    #     """
    #     print(f"\n{'=' * 25}\n--- Starting Full Atlas Analysis ---")
    #     print(f"Using Score: '{score_to_analyze}'")
    #     print(f"Using Distance Metric: '{scoring_system}'")
    #     print(f"{'=' * 25}\n")
    #
    #     # Clear any previous results to ensure a fresh run
    #     self.results = {}
    #
    #     # Ensure the results directory exists
    #     os.makedirs(results_path, exist_ok=True)
    #
    #     cancer_files = sorted(os.listdir(self.__cancer_csvs_path))
    #     if not cancer_files:
    #         print("[Warning] No cancer files found in the specified directory.")
    #         return {}
    #
    #     backup_filename = f"backup_results_{scoring_system}_{score_to_analyze}.pickle"
    #     backup_file_path = os.path.join(results_path, backup_filename)
    #
    #     for i, cancer_file in enumerate(cancer_files):
    #         cancer_file_path = os.path.join(self.__cancer_csvs_path, cancer_file)
    #
    #         # Skip directories or non-file items
    #         if not os.path.isfile(cancer_file_path):
    #             continue
    #
    #         print(f"\n{'#' * 15} Processing Cancer {i + 1}/{len(cancer_files)}: {cancer_file} {'#' * 15}")
    #
    #         # This method will perform the analysis for one cancer vs all pathways
    #         # and internally update self.results
    #         self.calculate_distances_single_cancer(
    #             cancer_file_path,
    #             score_to_analyze,
    #             scoring_system
    #         )
    #
    #         # --- Save a backup after each cancer is fully processed ---
    #         try:
    #             with open(backup_file_path, 'wb') as f:
    #                 pickle.dump(self.results, f, protocol=pickle.HIGHEST_PROTOCOL)
    #             print(f"\n[Backup] Progress for {i + 1} cancers saved to: {backup_file_path}")
    #         except Exception as e:
    #             print(f"\n[Error] Could not save backup file: {e}")
    #
    #     print(f"\n{'=' * 30}\n--- Full Atlas Analysis Completed ---\n{'=' * 30}")
    #     return self.results


"""
############ WORKING EXAMPLE: ############


# 1. Instantiate the class with the paths to your files.
#    Replace these paths with the actual locations of your files.
pathway_file = "./data/kegg/pathways/objects/hsa04010.pickle"
cancer_file = "/cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/data/cbio/cancers/acc_mutations.csv"
analyzer = CancerPathwayScoring(pathway_file, cancer_file)

# 2. Collect the background scores from all genes in the pathway.
print(f"Step 1: Collecting background scores for pathway '{pathway_file}'...")
background_scores = analyzer.get_pathway_scores_background()

# 3. Collect the cancer-specific scores for genes in that same pathway.
print(f"Step 2: Collecting cancer-specific scores from '{cancer_file}'...")
cancer_scores = analyzer.get_cancer_pathway_scores()

# 4. Proceed only if both datasets have scores to compare.
if background_scores and cancer_scores:
    # We will use the 'clinvar_reg_dis_ordered_prob' score for this example.
    score_to_analyze = 'clinvar_reg_dis_ordered_prob'

    # --- KL-Divergence Calculation ---
    print(f"\nStep 3a: Calculating KL-Divergence using '{score_to_analyze}' scores...")
    kl_divergence = analyzer.calculate_distance_gmm_kl_d(
        background_scores,
        cancer_scores,
        score_type=score_to_analyze
    )

    # --- Wasserstein Distance Calculation ---
    print(f"Step 3b: Calculating Wasserstein Distance using '{score_to_analyze}' scores...")
    wasserstein_dist = analyzer.calculate_distance_wasserstein(
        background_scores,
        cancer_scores,
        score_type=score_to_analyze
    )

    # 5. Print the final results.
    print("\n--- Analysis Complete ---")
    if kl_divergence is not None:
        print(f"The KL-Divergence is: {kl_divergence:.4f}")
        print("(Measures the information gain from switching from the background to the cancer model.)")
    else:
        print("Could not calculate KL-Divergence.")

    if wasserstein_dist is not None:
        print(f"The Wasserstein Distance is: {wasserstein_dist:.4f}")
        print("(Represents the 'work' required to transform the background distribution into the cancer distribution.)")
    else:
        print("Could not calculate Wasserstein Distance.")

    print(
        "\n(Higher values for both metrics indicate a greater difference between the cancer mutation profile and the background model.)")

else:
    print("\n--- Analysis Skipped ---")
    if not background_scores:
        print("Could not find any background scores in the pathway file.")
    if not cancer_scores:
        print("No mutations found in the cancer dataset for the genes in this specific pathway.")
"""