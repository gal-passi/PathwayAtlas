from definitions import *
import pandas as pd
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
import pickle




class BackgroundModel:
    def __init__(self):
        self.pathways_folder = KEGG_PATHWAY_OBJECTS_PATH

        # TODO change pssm scoring to michal`s pssm
        self.pssm = UNIFY_PSSM


    def collect_scores(self):
        """
        For each pathway, go through the genes in it and collect all mutation scoring:
        ["esm_min_max_naive" - min-max on score, "clinvar_reg_..." - from the 3 regressors],
        splitting them into the 12 mutation types.

        Returns
        -------
        dict:
            { pathway_name:
                { "A>C": {"esm_min_max_naive": [...],
                            "clinvar_reg_dis_ordered_prob": [...], "clinvar_reg_global_prob": [...]},
                  "A>G": {...}, ...
                }
            }
        """
        results = {}

        for file_name in os.listdir(self.pathways_folder)[0:1000]:     # TODO remove debug here [0:100]
            if not file_name.endswith(".pickle"):
                continue

            pathway_file = os.path.join(self.pathways_folder, file_name)

            try:
                with open(pathway_file, "rb") as f:
                    gene_to_csv = pickle.load(f)  # dict: gene_id -> csv path
            except Exception:
                continue  # skip unreadable

            pathway_scores, pathway_id = self.pathway_scores_collecting(file_name, gene_to_csv)

            if pathway_scores:
                results[pathway_id] = dict(pathway_scores)

        return results

    def pathway_scores_collecting(self, file_name, gene_to_csv):
        pathway_id = file_name.replace(".pickle", "")
        pathway_scores = defaultdict(lambda: {"esm_min_max_naive": [],
                                              "clinvar_reg_dis_ordered_prob": [],
                                              "clinvar_reg_global_prob": []})
        for gene_id, csv_path in gene_to_csv.items():
            if isinstance(csv_path, str) and os.path.exists(csv_path):
                try:
                    df = pd.read_csv(csv_path)
                except Exception:
                    continue

                # must have these cols
                if not {"Ref", "Alt"}.issubset(df.columns):
                    continue

                for _, row in df.iterrows():
                    ref, alt = str(row["Ref"]).upper(), str(row["Alt"]).upper()
                    mut_type = f"{ref}>{alt}"

                    for score_type in SCORE_TYPES:
                        try:
                            score_val = float(row[score_type])
                        except (ValueError, TypeError, KeyError):
                            score_val = None

                        pathway_scores[mut_type][score_type].append(score_val)

        return pathway_scores, pathway_id

    def create_joint_distribution(self, pathway_scores,
                                  score_type="clinvar_reg_dis_ordered_prob", bins=BINS_NUM, range=(0, 1)):
        """
        Calculate the joint distribution for a pathway using the pssm matrix as weights.

        This function combines the individual score distributions of each mutation type,
        weighting them by their frequency in self.pssm, to create a single
        mixed background distribution.

        Parameters
        ----------
        pathway_scores : dict
            A dictionary containing scores for each mutation type.
            Example: { "A>C": {"esm_min_max_naive": [...]}, ... }
        score_type : str, optional
            The type of score to use from pathway_scores, by default "esm_min_max_naive".
        bins : int, optional
            The number of bins to use for creating the histograms, by default BINS_NUM.
        range : tuple, optional
            The range of scores to consider for the histogram, by default (0, 1).

        Returns
        -------
        tuple:
            - np.ndarray: The values (heights) of the mixed distribution histogram.
            - np.ndarray: The bin edges for the histogram.
        """
        # Initialize an array to hold the final mixed distribution (histogram y-values).
        joint_distribution = np.zeros(bins)
        # To normalize correctly, we need to track the sum of weights for which we have data.
        total_weight_used = 0.0     # TODO the weight is the same all the time - the sum of the 12 relevant cells...

        # Pre-calculate bin edges to be used for all histograms.
        # This ensures all histograms are on the same x-axis.
        bin_edges = np.linspace(range[0], range[1], bins + 1)

        # Iterate over all possible mutation types and their a priori weights from the PSSM.
        for mut_type, weight in self.pssm.items():
            # Check if the pathway has data for this mutation type.
            if mut_type in pathway_scores and pathway_scores[mut_type][score_type]:

                # Get the scores and filter out invalid values (e.g., None, NaN).
                scores = pathway_scores[mut_type][score_type]
                valid_scores = [s for s in scores if s is not None and np.isfinite(s)]

                # If there are no valid scores for this mutation type, skip it.
                if not valid_scores:
                    continue

                # Calculate the histogram for the current mutation type's scores.
                # Using density=True normalizes the histogram so the area under it is 1,
                # turning it into a probability density estimate.
                hist, _ = np.histogram(valid_scores, bins=bin_edges, density=True)

                # Add the weighted distribution to our overall joint distribution.
                joint_distribution += hist * weight
                total_weight_used += weight

        # Normalize the final distribution by the sum of weights for the distributions
        # we actually included. This ensures that the final mixed distribution is a
        # proper probability distribution (integrates to 1).
        if total_weight_used > 0:
            joint_distribution /= total_weight_used     # TODO the weight is the same all the time - the sum of the 12 relevant cells...

        return joint_distribution, bin_edges




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

    def save_distribution(self, pathway_id, joint_distribution, bin_edges,
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

