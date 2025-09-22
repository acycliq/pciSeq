"""Statistical calculation utilities."""
import numpy as np
import pandas as pd
import numpy_groupies as npg
from typing import Tuple, Optional, Any, Union
import logging
import opt_einsum as oe
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

# Configure logging
ops_utils_logger = logging.getLogger(__name__)


def expected_covariance(scale_matrix, dof):
    """
        Calculate the expected covariance matrix from a scale matrix and degrees of freedom.

        Parameters
        ----------
        scale_matrix : np.ndarray
            Scale matrix of shape (C, d, d) where d must be 2 or 3
        dof : np.ndarray
            Degrees of freedom,shape (C,).
            Values will be automatically adjusted if below d + 2

        Returns
        -------
        np.ndarray
            Expected covariance matrix of same shape as input scale_matrix

        Raises
        ------
        ValueError
            If matrix dimensions are invalid or don't match
    """
    # Get the last two dimensions
    *_, d1, d2 = scale_matrix.shape

    # Check square
    if d1 != d2:
        raise ValueError(f"scale_matrix must be square, got shape {scale_matrix}")

    # Check dimension is 2 or 3
    if d1 not in (2, 3):
        raise ValueError(f"scale_matrix dimension must be 2 or 3, got {d1}")

    # Adjust degrees of freedom if needed, maybe I should drop a warning?
    min_dof = d1 + 1
    dof[dof <= min_dof] = min_dof + 1

    return scale_matrix / (dof[:, None, None] - d1 - 1)


def negative_binomial_loglikelihood(x: np.ndarray, r: float, q: np.ndarray) -> np.ndarray:
    """Calculate the Negative Binomial log-likelihood for given parameters.

    The Negative Binomial distribution models the number of failures (x) before
    observing the r-th success, with failure probability q. The PMF is:
        P(X = x) = C(x + r - 1, x) * q^x * (1 - q)^r

    Here we compute only the terms that depend on q and r:
        log-likelihood = x * log(q) + r * log(1 - q)

    Args:
        x: Array of observed failure counts (non-negative floats).
        r: Number of successes until stopping (dispersion parameter, positive).
        q: Array of failure probabilities (each between 0 and 1).

    Returns:
        Array of log-likelihood values, broadcast over x and q.

    Raises:
        ValueError: If any q is outside (0, 1) or if x has negative values.
    """
    try:
        x = x[:, :, None]  # Add dimension for broadcasting

        # Compute the log-likelihood of seeing x failures before the r-th success,
        # if the failure probability is q.
        # In our context, x is the cell gene counts, q is derived from the single cell data
        # count data and r is a hyperparameter (set by default = 2.0).
        # Scipy's nbinom object has logpmf(k, n, p) where p is the prob of success, ie p = 1-q
        # and k, n is what is denoted here by x, r respectively. Also logpmf includes the
        # combinatorial factor. Finally logpmf will drop an exception if the counts k are not
        # integers
        log_likelihood = x * np.log(q) + r * np.log(1 - q)

        return log_likelihood

    except Exception as e:
        ops_utils_logger.error(f"Error calculating negative binomial log-likelihood: {str(e)}")
        raise ValueError("Failed to compute log-likelihood. Check input dimensions and values.")


def compute_gene_loglikelihood_matrix(obj) -> np.ndarray:
    """
    Compute the full gene log-likelihood contribution matrix for all cells and cell types.
    This function performs the core computation shared between cell_to_cellType and
    calculate_genes_log_likelihood_contr, eliminating code duplication and improving performance.
    Args:
        obj: VarBayes object containing the following attributes:
            - scaled_exp: A delayed or computed array of scaled expression values (shape: nC x nG x nK)
            - genes.eta_bar: Gene efficiency parameters (shape: nG)
            - config['SpotReg']: Regularization parameter for spot-level noise
            - config['rSpot']: Dispersion parameter for the negative binomial distribution
            - cells.geneCount: Observed gene counts for all cells (shape: nC x nG)
    Returns:
        np.ndarray: Log-likelihood contributions matrix of shape (nC, nG, nK)
                   where element [c,g,k] is the log-likelihood contribution of
                   gene g in cell c for cell type k
    """
    # Compute scaled expression (expensive operation done once)
    scaled_means = obj.scaled_exp.compute()

    # Calculate scaled expression adjusted by gene efficiency and regularization
    ScaledExp = np.einsum('cgk,cg,g->cgk',
                          scaled_means,
                          obj.cells.plane_adj.values,
                          obj.genes.eta_bar) + obj.config['SpotReg']

    # Calculate negative binomial probabilities
    pNegBin = ScaledExp / (obj.config['rSpot'] + ScaledExp)

    # Get gene counts for all cells
    cgc = obj.cells.geneCount

    # Calculate log-likelihood contributions for all cells
    contr = negative_binomial_loglikelihood(cgc, obj.config['rSpot'], pNegBin)

    return contr



def calculate_genes_log_likelihood_contr(obj, label: int) -> Tuple[DataFrame, Series, DataFrame]:
    """
    Calculate the log-likelihood contributions, gene counts, and scaled expression values
    for a specific cell.

    This function computes:
        1. The genes' log-likelihood contributions (`contr`) for the specified cell under a
           negative binomial distribution.
        2. The gene counts (`cgc`) for the specified cell.
        3. The scaled expression values (`scaled_means`) for the specified cell.

    Args:
        obj: An object containing the following attributes:
            - scaled_exp: A delayed or computed array of scaled expression values (shape: nC x nG x nK).
            - genes.eta_bar: Gene efficiency parameters (shape: nG).
            - config['SpotReg']: Regularization parameter for spot-level noise.
            - config['rSpot']: Dispersion parameter for the negative binomial distribution.
            - cells.geneCount: Observed gene counts for all cells (shape: nC x nG).
        label (int): The index of the cell for which to compute the values.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
            - contr: The log-likelihood contributions for the specified cell (shape: nG x nK).
            - cgc: The gene counts for the specified cell (shape: nG).
            - scaled_means: The scaled expression values for the specified cell (shape: nG x nK).
    """
    # If original labels have been renumbered find the label it's been mapped to.
    if obj.config['label_map']:
        label = obj.config['label_map'][label]

    # Get the full log-likelihood matrix using shared computation
    contr = compute_gene_loglikelihood_matrix(obj)

    # Get scaled expression and gene counts
    scaled_means = obj.scaled_exp.compute()

    # Get gene counts for all cells
    cgc = obj.cells.geneCount

    # Return values for the specified cell
    contr_df = pd.DataFrame(contr[label], columns=obj.cells.class_names).set_index(obj.genes.gene_panel)
    gene_counts = pd.Series(cgc[label], index=obj.genes.gene_panel)
    scaled_means_df = pd.DataFrame(scaled_means[label], columns=obj.cells.class_names).set_index(obj.genes.gene_panel)
    return contr_df, gene_counts, scaled_means_df


# def plot_loglik_contr(df):
#     """
#     Create a scatter plot of the first column vs the second column in a DataFrame,
#     with tooltips from the index, and add a diagonal line (y = x).
#
#     Args:
#         df (pd.DataFrame): The DataFrame containing the data.
#     """
#     # Ensure the DataFrame has at least two columns
#     if len(df.columns) < 2:
#         raise ValueError("The DataFrame must have at least two columns.")
#
#     # Reset the index to include it as a column for tooltips
#     df = df.reset_index()
#
#     # Get the names of the first and second columns
#     x_col = df.columns[1]  # First column (after resetting the index)
#     y_col = df.columns[2]   # Second column (after resetting the index)
#
#     # Create the scatter plot with tooltips
#     fig = px.scatter(
#         df,
#         x=x_col,
#         y=y_col,
#         hover_data=['index'],  # Include the index as a tooltip
#         title=f"Scatter Plot: {x_col} vs {y_col}"
#     )
#
#     # Add a diagonal line (y = x)
#     min_val = min(df[x_col].min(), df[y_col].min())  # Minimum value across both axes
#     max_val = max(df[x_col].max(), df[y_col].max())  # Maximum value across both axes
#
#     diagonal_line = go.Scatter(
#         x=[min_val, max_val],  # X values for the line (y = x)
#         y=[min_val, max_val],  # Y values for the line (y = x)
#         mode='lines',  # Draw a line
#         name='Diagonal Line (y = x)',  # Label for the line
#         line=dict(color='red', dash='dash')  # Customize line color and style
#     )
#
#     # Add the diagonal line to the figure
#     fig.add_trace(diagonal_line)
#
#     # Show the plot
#     fig.show()


# def visualize_fit(gene_counts, scaled_means):
#     """
#     Visualize the fit between gene_counts and scaled_means using Plotly.
#
#     Args:
#         gene_counts (pd.Series): Observed gene counts for a cell.
#         scaled_means_df (pd.DataFrame): Scaled expected gene expression values for the cell.
#     """
#     # Ensure gene_counts and scaled_means_df have the same index (gene names)
#     if not gene_counts.index.equals(scaled_means.index):
#         raise ValueError("gene_counts and scaled_means_df must have the same index.")
#
#     for column in scaled_means.columns:
#         # Create a scatter plot
#         fig = go.Figure()
#
#         # Add scatter plot: gene_counts vs. scaled_means
#         scatter_trace = go.Scatter(
#             x=scaled_means[column],
#             y=gene_counts,
#             mode='markers',
#             marker=dict(opacity=0.6),
#             text=gene_counts.index,  # Tooltip: gene names
#             name='Scatter Plot'
#         )
#         fig.add_trace(scatter_trace)
#
#         # Add a true diagonal line (y = x)
#         min_val = min(scaled_means[column].min(), gene_counts.min())  # Minimum value across both axes
#         max_val = max(scaled_means[column].max(), gene_counts.max())  # Maximum value across both axes
#
#         diagonal_line = go.Scatter(
#             x=[min_val, max_val],  # X values for the line (y = x)
#             y=[min_val, max_val],  # Y values for the line (y = x)
#             mode='lines',
#             line=dict(color='red', dash='dash'),
#             name='y = x'
#         )
#         fig.add_trace(diagonal_line)
#
#         # Update layout
#         fig.update_layout(
#             title=f'Gene Counts vs. Scaled Means ({column})',
#             xaxis_title=f'Scaled Means ({column})',
#             yaxis_title='Gene Counts',
#             showlegend=True
#         )
#
#         # Calculate correlation
#         correlation = gene_counts.corr(scaled_means[column])
#
#         # Calculate residuals and their sum
#         residuals = gene_counts - scaled_means[column]
#         sum_residuals = residuals.sum()
#
#         # Print correlation and sum of residuals
#         print(f"Correlation between gene_counts and {column}: {correlation:.3f}")
#         print(f"Sum of residuals for {column}: {sum_residuals:.3f}")
#
#         # Show the plot
#         fig.show()


def check_cell(obj, label, user_class, top_n=10, show_plot=True):
    """
    Compare gene expression likelihoods between two classes for a specific cell.

    Parameters:
        label (int): The cell number to analyze.
        user_class (str): The user-specified class to compare against.
        top_n (int): Number of top and bottom genes to retrieve (default: 10).

    Returns:
        pd.DataFrame: A DataFrame containing mean expression values and gene counts for the top and bottom genes.
    """

    # If original labels have been renumbered find the label it's been mapped to.
    if obj.config['label_map']:
        pciSeq_label = obj.config['label_map'][label]
    else:
        pciSeq_label = label

    # Step 1: Calculate gene log-likelihood contributions
    contr_df, gene_counts, _ = obj.calculate_genes_log_likelihood_contr(label)

    # Step 2: Get the cell's class from cellData
    pciSeq_class = obj.cells.class_names[obj.cells.classProb[pciSeq_label].argmax()]

    # Step 3: Check if classes exist in contr_df
    if pciSeq_class not in contr_df.columns or user_class not in contr_df.columns:
        raise ValueError(f"One or both classes ({pciSeq_class}, {user_class}) not found in contr_df.")

    # Step 4: Calculate differences and get top/bottom genes
    my_contr_df = contr_df[[pciSeq_class, user_class]].copy()
    my_contr_df['diff'] = my_contr_df[pciSeq_class] - my_contr_df[user_class]

    top_genes = my_contr_df.nlargest(top_n, 'diff').index.values
    bottom_genes = my_contr_df.nsmallest(top_n, 'diff').index.values

    # Step 5: Combine top and bottom genes
    selected_genes = np.append(top_genes, bottom_genes)

    # Step 6: Retrieve mean expression and gene counts
    # gene_expression_data = obj.single_cell.mean_expression.loc[selected_genes, [pciSeq_class, user_class]]
    # gene_expression_data = gene_expression_data.merge(
    #     gene_counts[selected_genes].rename('Cell Gene Counts'),
    #     left_index=True,
    #     right_index=True
    # )

    # Step 6: Retrieve mean expression and gene counts
    gene_expression_data = pd.DataFrame(
        obj.cells.mean_gene_reads_per_class(),
        columns=obj.cells.class_names
    ).set_index(obj.genes.gene_panel)

    # Filter rows and columns
    gene_expression_data = gene_expression_data.loc[selected_genes, [pciSeq_class, user_class]]

    # Merge with gene_counts
    gene_expression_data = gene_expression_data.merge(
        gene_counts[selected_genes].rename('Cell Gene Counts'),
        left_index=True,
        right_index=True
    )

    # Add the MultiIndex header
    new_columns = pd.MultiIndex.from_tuples(
        [('pciSeq class avg counts', pciSeq_class),
         ('pciSeq class avg counts', user_class),
         (f'Cell: {label}', 'Gene Counts')]  # Leave the third column header as is
    )
    gene_expression_data.columns = new_columns

    if show_plot:
        # Step 7: Plot top and bottom genes as subplots
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Calculate the sum of the top n contributions
        top_contribution_sum = my_contr_df.loc[top_genes, 'diff'].sum()
        bottom_contribution_sum = my_contr_df.loc[bottom_genes, 'diff'].sum()

        # Plot top genes
        my_contr_df.loc[top_genes, 'diff'].plot.bar(ax=axes[0], color='skyblue',
                                                    title=f'Cell: {label} - Top {top_n} contr for class: {pciSeq_class} (Sum: {top_contribution_sum:.2f})')
        axes[0].set_ylabel('Log-Likelihood Difference')
        axes[0].set_xlabel('Genes')

        # Plot bottom genes
        my_contr_df.loc[bottom_genes, 'diff'].plot.bar(ax=axes[1], color='lightcoral',
                                                       title=f'Cell: {label} - Top {top_n} contr for class: {user_class} (Sum: {bottom_contribution_sum:.2f})')
        axes[1].set_ylabel('Log-Likelihood Difference')
        axes[1].set_xlabel('Genes')

        plt.tight_layout()
        plt.show()

    return gene_expression_data, fig if show_plot else None


def read_tsv(filepath):
    """
    Convenience function to read the tsv files generated by pciSeq
    """
    data = pd.read_csv(filepath, sep='\t')
    data = data.map(
        lambda x: eval(x) if isinstance(x, str) and x.strip().startswith(('{', '[', '(')) else x)
    return data


# def softmax(X: np.ndarray, theta: float = 1.0, axis: Optional[int] = None) -> np.ndarray:
#     """Compute the softmax of each element along an axis of X.
#
#     Args:
#         X: Input array (should be floats)
#         theta: Multiplier prior to exponentiation (default: 1.0)
#         axis: Axis to compute values along (default: first non-singleton axis)
#
#     Returns:
#         Array same size as X, normalized along the specified axis
#
#     Notes:
#         From https://nolanbconaway.github.io/blog/2017/softmax-numpy
#     """
#     # Make X at least 2d
#     y = np.atleast_2d(X)
#
#     # Find axis if not specified
#     if axis is None:
#         axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)
#
#     # Multiply y against the theta parameter
#     y = y * float(theta)
#
#     # Subtract the max for numerical stability
#     y = y - np.expand_dims(np.max(y, axis=axis), axis)
#
#     # Exponentiate y
#     y = np.exp(y)
#
#     # Take the sum along the specified axis
#     ax_sum = np.expand_dims(np.sum(y, axis=axis), axis)
#
#     # Finally: divide elementwise
#     p = y / ax_sum
#
#     # Flatten if X was 1D
#     if len(X.shape) == 1:
#         p = p.flatten()
#
#     return p


def has_converged(
        spots: Any,
        p0: Optional[np.ndarray],
        tol: float
) -> Tuple[bool, float]:
    """Check if probability assignments have converged.

    Args:
        spots: Spot data object containing parent_cell_prob
        p0: Previous probability matrix (None for first iteration)
        tol: Convergence tolerance threshold

    Returns:
        Tuple containing:
            - bool: True if converged, False otherwise
            - float: Maximum absolute difference between iterations

    Raises:
        Exception: If convergence check fails
    """
    p1 = spots.parent_cell_prob
    if p0 is None:
        p0 = np.zeros_like(p1)

    try:
        delta = np.max(np.abs(p1 - p0))
        converged = (delta < tol)
        return converged, delta
    except Exception as e:
        ops_utils_logger.error(f"Convergence check failed: {str(e)}")
        raise


def scaled_exp(cell_area_factor: np.ndarray,
               sc_mean_expressions: np.ndarray) -> np.ndarray:
    """Calculate scaled expression values.

    Args:
        cell_area_factor: Cell area scaling factors
        sc_mean_expressions: Single cell mean expression values

    Returns:
        Scaled expression array
    """
    subscripts = 'c,gk->cgk'
    operands = [cell_area_factor, sc_mean_expressions]

    return oe.contract(subscripts, *operands, optimize='optimal')


def empirical_mean(spots, cells):

    # get the total gene counts per cell
    N_c = cells.total_counts

    xyz_spots = spots.xyz_coords
    prob = spots.parent_cell_prob
    n = cells.config['nNeighbors'] + 1

    # multiply the x coord of the spots by the cell prob
    a = np.tile(xyz_spots[:, 0], (n, 1)).T * prob

    # multiply the y coord of the spots by the cell prob
    b = np.tile(xyz_spots[:, 1], (n, 1)).T * prob

    # multiply the z coord of the spots by the cell prob
    c = np.tile(xyz_spots[:, 2], (n, 1)).T * prob

    # aggregated x and y coordinate
    idx = spots.parent_cell_id
    x_agg = npg.aggregate(idx.ravel(), a.ravel(), size=len(N_c))
    y_agg = npg.aggregate(idx.ravel(), b.ravel(), size=len(N_c))
    z_agg = npg.aggregate(idx.ravel(), c.ravel(), size=len(N_c))

    # get the estimated cell centers
    x_bar = np.nan * np.ones(N_c.shape)
    y_bar = np.nan * np.ones(N_c.shape)
    z_bar = np.nan * np.ones(N_c.shape)

    x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
    y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]
    z_bar[N_c > 0] = z_agg[N_c > 0] / N_c[N_c > 0]

    # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
    xyz_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T, z_bar.T)))

    # if you have a value for the estimated centroid use that, otherwise
    # use the initial (starting values) centroids
    ini_cent = cells.ini_centroids()
    xyz_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y'], ini_cent['z']])))

    # # sanity check. NaNs or Infs should appear together
    # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
    # use the fitted centroids where possible otherwise use the initial ones
    xyz_bar[np.isfinite(x_bar)] = xyz_bar_fitted[np.isfinite(x_bar)]
    return pd.DataFrame(xyz_bar, columns=['x', 'y', 'z'], dtype=np.float32)


def gene_density(spots, config) -> pd.DataFrame:
    """
    Calculate gene density adjustment values across imaging planes.

    PROBLEM: In our pipeline gene detection efficiency varies
    across z-planes. Some planes may have very few reads for certain genes
    while others have abundant reads for the same genes, creating detection bias
    that can lead to incorrect cell type classification.

    SOLUTION: Compute plane-specific adjustment values for each gene by
    comparing the gene's expression in each plane to its average expression across
    all planes. The resulting density matrix can be used to adjust single-cell
    expression data to account for plane-specific detection efficiency variations.

    The normalization process:
    1. Counts spots per (plane, gene) combination
    2. Calculates mean counts per gene across planes (excluding zeros)
    3. Normalizes each plane's counts by the gene's mean
    4. Sets density to 1.0 for genes absent in a plane (no adjustment)

    A value of 1.o means that no adjustment is applied to the single-cell data.
    The value post-adjustment is the same as the original value

    COMMENT to myself:Maybe I should also introduce a regularisation parameter too
    """

    data = spots.data.assign(gene_id=spots.gene_id)

    # Count spots per plane/gene_name and pivot
    counts = data.groupby(["plane_id", "gene_name"]).size().unstack(fill_value=0)

    # Ensure all planes are represented
    all_planes = np.arange(config['img_dim']['n_planes'])
    counts = counts.reindex(index=all_planes, columns=sorted(counts.columns), fill_value=0)

    # Calculate means over non-zero values only
    gene_means = counts.replace(0, np.nan).mean(axis=0)

    # Normalize by gene means (density calculation).
    # if an element is zero, set it to 1. Effectively that means that if the gene is not present on that plane,
    # then dont make any adjustment to the single cell data
    density = counts.div(gene_means, axis=1).fillna(0).astype(np.float32)
    density[density == 0] = 1

    # print("REMOVE THIS - REMOVE THIS")
    # density = pd.DataFrame(np.ones(density.shape)) # REMOVE THIS - REMOVE THIS
    return density # num_planes x num_genes


