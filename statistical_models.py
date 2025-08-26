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
        self.pathways_and_modules_folder = KEGG_PATHWAY_OBJECTS_PATH

        # TODO find pssm scoring
        # TODO move to definitions
        self.pssm = {
            "A>G": 1 / 12, "A>C": 1 / 12, "A>T": 1 / 12,
            "G>A": 1 / 12, "G>C": 1 / 12, "G>T": 1 / 12,
            "T>A": 1 / 12, "T>C": 1 / 12, "T>G": 1 / 12,
            "C>A": 1 / 12, "C>G": 1 / 12, "C>T": 1 / 12
        }


    def collect_scores(self):
        """
        For each pathway, go through the genes in it and collect all mutation scores:
        ["normalized_score", "esm_disorder_scoring"],
        splitting them into the 12 mutation types.

        Returns
        -------
        dict:
            { pathway_name:
                { "A>C": {"normalized_score": [...], "esm_disorder_scoring": [...]},
                  "A>G": {...}, ...
                }
            }
        """
        results = {}

        for file_name in os.listdir(self.pathways_and_modules_folder):
            if not file_name.endswith(".pickle"):
                continue

            pathway_file = os.path.join(self.pathways_and_modules_folder, file_name)

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
        pathway_scores = defaultdict(lambda: {"normalized_score": [], "esm_disorder_scoring": []})

        for gene_id, csv_path in gene_to_csv.items():
            if not os.path.exists(csv_path):
                continue

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

                try:
                    ns = float(row["normalized_score"])
                except (ValueError, TypeError, KeyError):
                    ns = None

                try:
                    es = float(row["esm_disorder_scoring"])
                except (ValueError, TypeError, KeyError):
                    es = None

                pathway_scores[mut_type]["normalized_score"].append(ns)
                pathway_scores[mut_type]["esm_disorder_scoring"].append(es)

        return pathway_scores, pathway_id

    def create_joint_distribution(self, pathway_scores, score_type="normalized_score", bins=100, range=(0, 1)):
        """
        Calculate the joint distribution for a pathway using the pssm matrix as weights.

        This function combines the individual score distributions of each mutation type,
        weighting them by their frequency in self.pssm, to create a single
        mixed background distribution.

        Parameters
        ----------
        pathway_scores : dict
            A dictionary containing scores for each mutation type.
            Example: { "A>C": {"normalized_score": [...]}, ... }
        score_type : str, optional          [normalized_score or esm_disorder_scoring]
            The type of score to use from pathway_scores, by default "normalized_score".
        bins : int, optional
            The number of bins to use for creating the histograms, by default 100.
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
        total_weight_used = 0.0

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
            joint_distribution /= total_weight_used

        return joint_distribution, bin_edges

    def GMM_the_distribution(self, joint_distribution, bin_edges, max_components=10, n_samples=100000,
                             random_state=None):
        """
        Fit a Gaussian Mixture Model (GMM) to a distribution defined by a histogram.

        This method first reconstructs a dataset from the histogram data, then
        iteratively fits GMMs with different numbers of components. It uses the
        Bayesian Information Criterion (BIC) to select the optimal number of
        components and returns the best-fitted model. [1, 5, 6]

        Parameters
        ----------
        joint_distribution : np.ndarray
            The values (heights) of the histogram representing the distribution.
        bin_edges : np.ndarray
            The edges of the histogram bins.
        max_components : int, optional
            The maximum number of Gaussian components to test, by default 10.
        n_samples : int, optional
            The number of samples to generate for reconstructing the dataset from
            the histogram, by default 100000.
        random_state : int, optional
            A random state for reproducibility, by default None.

        Returns
        -------
        sklearn.mixture.GaussianMixture
            The best fitted Gaussian Mixture Model object, selected based on the
            lowest BIC score. [3] Returns None if the input distribution is empty.
        """
        # 1. Reconstruct the dataset from the histogram
        # Calculate the center of each bin
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Calculate the number of samples to draw from each bin, proportional to its height
        # This creates a representative dataset from the probability distribution
        counts = (joint_distribution * n_samples).astype(int)

        # If there's no data to model, return None
        if np.sum(counts) == 0:
            return None

        # Repeat each bin center 'count' times to generate the data
        reconstructed_data = np.repeat(bin_centers, counts).reshape(-1, 1)

        # 2. Find the optimal number of components using BIC
        lowest_bic = np.inf
        best_gmm = None

        # Define the range of component numbers to test
        n_components_range = range(1, max_components + 1)

        for n_components in n_components_range:
            # Fit a Gaussian Mixture Model
            gmm = GaussianMixture(n_components=n_components,
                                  random_state=random_state,
                                  covariance_type='full')
            gmm.fit(reconstructed_data)

            # Calculate the BIC for the current model. [4]
            bic = gmm.bic(reconstructed_data)

            # If the current model has a lower BIC, it's a better fit
            if bic < lowest_bic:
                lowest_bic = bic
                best_gmm = gmm

        # 3. Return the best model found
        return best_gmm

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

        return ax

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

    # TODO add KDE smoothing option? It might be better in our case
