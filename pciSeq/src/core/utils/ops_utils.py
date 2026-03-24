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
from scipy.special import psi, softmax

# Configure logging
logger = logging.getLogger(__name__)


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
        logger.error(f"Error calculating negative binomial log-likelihood: {str(e)}")
        raise ValueError("Failed to compute log-likelihood. Check input dimensions and values.")


def compute_gene_loglikelihood_matrix(obj) -> np.ndarray:
    """
    Compute the full gene log-likelihood contribution matrix for all cells and cell types.

    This function performs the core computation shared between cell_to_cellType and
    calculate_genes_log_likelihood_contr, eliminating code duplication and improving performance.

    Args:
        obj: VarBayes object containing the following attributes:
            - scaled_exp: A delayed or computed array of scaled expression values (shape: nC x nG x nK)
            - genes.eta_bar: Gene efficiency (shape: nG)
            - cells.theta_bar: Cell inefficiency (shape: nC)
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
    ScaledExp = np.einsum('cgk,g,ck->cgk', scaled_means, obj.genes.eta_bar, obj.cells.theta_bar) + obj.config['SpotReg']

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
    new_columns = pd.MultiIndex.from_tuples([
        (f'Cells typed as {pciSeq_class}', 'mean counts'),
        (f'Cells typed as {user_class}', 'mean counts'),
        (f'This cell: ({label})', 'counts')
    ])
    gene_expression_data.columns = new_columns

    # Step 7: Compute the prior and MRF terms for the two classes
    class_names = list(obj.cells.class_names)
    pciSeq_idx = class_names.index(pciSeq_class)
    user_idx = class_names.index(user_class)

    log_prior = obj.cellTypes.log_prior
    mrf = obj.cells.calc_mrf()

    gene_loglik_pciSeq = my_contr_df[pciSeq_class].sum()
    gene_loglik_user = my_contr_df[user_class].sum()
    log_prior_pciSeq = log_prior[pciSeq_idx]
    log_prior_user = log_prior[user_idx]
    mrf_pciSeq = mrf[pciSeq_label, pciSeq_idx]
    mrf_user = mrf[pciSeq_label, user_idx]

    # Log-posterior for the two classes
    log_post_pciSeq = gene_loglik_pciSeq + log_prior_pciSeq + mrf_pciSeq
    log_post_user = gene_loglik_user + log_prior_user + mrf_user

    # Posterior probabilities (softmax over just these two classes)
    log_posts = np.array([log_post_pciSeq, log_post_user])
    posterior_probs = softmax(log_posts)

    if show_plot:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # --- Top row: gene-level log-likelihood differences (unchanged) ---
        top_contribution_sum = my_contr_df.loc[top_genes, 'diff'].sum()
        bottom_contribution_sum = my_contr_df.loc[bottom_genes, 'diff'].sum()

        my_contr_df.loc[top_genes, 'diff'].plot.bar(ax=axes[0, 0], color='skyblue',
                                                    title=f'Cell: {label} - Top {top_n} contr for class: {pciSeq_class} (Sum: {top_contribution_sum:.2f})')
        axes[0, 0].set_ylabel('Log-Likelihood Difference')
        axes[0, 0].set_xlabel('Genes')

        my_contr_df.loc[bottom_genes, 'diff'].plot.bar(ax=axes[0, 1], color='lightcoral',
                                                       title=f'Cell: {label} - Top {top_n} contr for class: {user_class} (Sum: {bottom_contribution_sum:.2f})')
        axes[0, 1].set_ylabel('Log-Likelihood Difference')
        axes[0, 1].set_xlabel('Genes')

        # --- Bottom-left: grouped bar chart of log-posterior components ---
        x = np.arange(3)
        width = 0.35
        vals_pciSeq = [gene_loglik_pciSeq, log_prior_pciSeq, mrf_pciSeq]
        vals_user = [gene_loglik_user, log_prior_user, mrf_user]

        axes[1, 0].bar(x - width/2, vals_pciSeq, width, label=pciSeq_class, color='skyblue')
        axes[1, 0].bar(x + width/2, vals_user, width, label=user_class, color='lightcoral')
        axes[1, 0].set_xticks(x)
        axes[1, 0].set_xticklabels(['Gene LogLik', 'Log Prior', 'MRF'])
        axes[1, 0].set_ylabel('Log-scale value')
        axes[1, 0].set_title(f'Cell: {label} - Log-posterior components')
        axes[1, 0].legend()
        axes[1, 0].axhline(y=0, color='grey', linestyle='--', linewidth=0.5)

        # --- Bottom-right: posterior probabilities ---
        axes[1, 1].bar([pciSeq_class, user_class],
                       [posterior_probs[0] * 100, posterior_probs[1] * 100],
                       color=['skyblue', 'lightcoral'])
        axes[1, 1].set_ylabel('Posterior Probability (%)')
        axes[1, 1].set_title(f'Cell: {label} - Posterior probabilities')

        plt.tight_layout()
        plt.show()

    return gene_expression_data, my_contr_df, fig if show_plot else None


def cell_typing_breakdown(obj, label, weights=None, show_plot=True):
    """
    Follow cell-typing step-by-step for a given cell and assuming spot assignment is known

    Parameters:
        obj: The VarBayes object
        label (int): The cell label to analyze
        weights: Optional override for initial Dirichlet alpha.
            - dict: Same semantics as config['cell_type_weights']
              {'default': value_1, 'Class_1': value_2, ..., 'Class_n': value_n}.
              Unknown class keys are ignored with a warning. Values map by name
              to the order in obj.cellTypes.names.
            - 1D array-like: Explicit alpha vector of length K matching
              obj.cellTypes.names order. When using an array, you must include
              the entry for the 'Zero' class yourself.
        show_plot (bool): Whether to display plots (default: True)

    Returns:
        dict: Contains all intermediate values and final probabilities
    """

    # Get configuration
    prior_mode = obj.config.get('cell_type_prior', 'uniform')

    if prior_mode != 'weighted':
        logger.warning(
            f"Function available only for 'weighted' cell type prior mode."
        )
        return dict()

    # Step 1: Get initial alpha (from config weights) or override
    def _build_alpha_from_dict(dct, names):
        """Build alpha vector following the same logic as cell_type_weights.

        - Start from default=1 (or provided)
        - Override per-class entries when present
        - Ignore unknown keys with a warning
        """
        default_val = dct.get('default', 1)
        # Initialize with defaults
        vals = {name: default_val for name in names}
        # Apply overrides
        for key, val in dct.items():
            if key == 'default':
                continue
            if key not in names:
                logger.warning(
                    f"Cell type '{key}' in weights dict not found in cell type names. Ignoring.")
                continue
            vals[key] = val
        # Handle Zero if not explicitly provided
        # if 'Zero' not in dct:
        #     non_zero_names = [n for n in names if n != 'Zero']
        #     vals['Zero'] = float(np.sum([vals[n] for n in non_zero_names]))
        # Return in the exact order of names
        return np.array([float(vals[n]) for n in names], dtype=float)

    names = obj.cellTypes.names
    nK = obj.nK # number of classes (aka cell types) including 'Zero'


    if weights is not None:
        if isinstance(weights, dict):
            ini_alpha = _build_alpha_from_dict(weights, names)
            alpha_source = 'override'
        else:
            ini_alpha = np.asarray(weights, dtype=float)
            if ini_alpha.shape != (nK,):
                raise ValueError(f"weights must have shape ({nK},), got {ini_alpha.shape}")
            alpha_source = 'override'
    else:
        ini_alpha = obj.cellTypes.ini_alpha()
        alpha_source = 'default'

    # Step 2: Get observed class sizes (zeta). This is basically the number of cells in each class.
    zeta = obj.cells.classProb.sum(axis=0)

    # WARNING: DUPLICATED CODE. Steps 3 and 4 below are already in cellClass.
    # If I change something in CellClass, I need to change it here too.
    # It is OK for now, but If we develop cellClass any further this will be a problem.

    # Step 3: Updated alpha (what dalpha_upd does)
    updated_alpha = zeta + ini_alpha

    # Step 4: Compute log_prior from updated alpha
    if obj.single_cell.isMissing or prior_mode == 'weighted':
        log_prior = psi(updated_alpha) - psi(updated_alpha.sum())
    else:
        prior = updated_alpha / updated_alpha.sum()
        log_prior = np.log(prior)

    # Step 5: Get gene log-likelihood for this cell
    contr_df, _, _ = calculate_genes_log_likelihood_contr(obj, label)
    gene_loglik = contr_df.sum(axis=0).values  # Sum over genes

    # Step 6: Compute log posterior
    log_posterior = gene_loglik + log_prior

    # Step 7: Apply softmax to get final probabilities
    posterior_probs = softmax(log_posterior)

    # Store results
    out = {
        'label': label,
        'cell_type_names': obj.cellTypes.names,
        'ini_alpha': ini_alpha,
        'zeta': zeta,
        'updated_alpha': updated_alpha,
        'log_prior': log_prior,
        'gene_loglik': gene_loglik,
        'log_posterior': log_posterior,
        'posterior_probs': posterior_probs,
        'prior_mode': prior_mode,
        'alpha_source': alpha_source,
        'predicted_class': obj.cellTypes.names[np.argmax(posterior_probs)],
        'predicted_prob': np.max(posterior_probs)
    }

    if show_plot:
        _plot_classification_steps(out)

    return out


def _plot_classification_steps(data):
    """Helper function to plot the classification trace."""

    from plotly.subplots import make_subplots

    cell_type_names = data['cell_type_names']
    n_types = len(cell_type_names)

    # Compute cell class prior (softmax of log_prior)
    cell_class_prior = softmax(data['log_prior'])

    # Create subplots: 4 rows x 2 columns (leave last slot empty)
    fig = make_subplots(
        rows=4, cols=2,
        subplot_titles=(
            '<b>Step 1: Initial Alpha</b><br><sub>(from config weights)</sub>',
            '<b>Step 2: Updated Alpha</b><br><sub>(ini_alpha + zeta)</sub>',
            '<b>Step 3: Cell Class Log Prior</b><br><sub>(from updated alpha)</sub>',
            '<b>Step 4: Cell Class Prior</b><br><sub>(softmax of log prior)</sub>',
            '<b>Step 5: Cell Class Log-Likelihood</b><br><sub>(from gene expression data)</sub>',
            '<b>Step 6: Cell Class Log Posterior</b><br><sub>(log-likelihood + log prior)</sub>',
            '<b>Step 7: Cell Class Posterior</b><br><sub>(softmax of log posterior)</sub>',
            ''  # Empty placeholder
        ),
        # Reduce spacing to make each subplot taller (same overall size)
        vertical_spacing=0.08,
        # Slightly increase space between left and right columns
        horizontal_spacing=0.12
    )

    # Color scheme
    colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    bar_colors = [colors[i % len(colors)] for i in range(n_types)]

    # Plot 1: Initial alpha (Row 1, Col 1)
    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=data['ini_alpha'],
        marker_color=bar_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>ini_alpha: %{y:.2f}<extra></extra>'
    ), row=1, col=1)

    # Plot 2: Updated alpha (Row 1, Col 2)
    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=data['updated_alpha'],
        marker_color=bar_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>updated_alpha: %{y:.2f}<extra></extra>'
    ), row=1, col=2)

    # Plot 3: Cell Class Log Prior (Row 2, Col 1)
    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=data['log_prior'],
        marker_color=bar_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>log_prior: %{y:.3f}<extra></extra>'
    ), row=2, col=1)

    # Plot 4: Cell Class Prior (Row 2, Col 2)
    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=cell_class_prior * 100,
        marker_color=bar_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>prior: %{y:.2f}%<extra></extra>'
    ), row=2, col=2)

    # Plot 5: Cell Class Log-Likelihood (Row 3, Col 1)
    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=data['gene_loglik'],
        marker_color=bar_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>log_likelihood: %{y:.1f}<extra></extra>'
    ), row=3, col=1)

    # Plot 6: Cell Class Log Posterior (Row 3, Col 2)
    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=data['log_posterior'],
        marker_color=bar_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>log_posterior: %{y:.1f}<extra></extra>'
    ), row=3, col=2)

    # Plot 7: Cell Class Posterior (Row 4, Col 1) with highlight for winner
    max_idx = np.argmax(data['posterior_probs'])
    final_colors = [colors[i % len(colors)] if i != max_idx else '#e74c3c'
                   for i in range(n_types)]

    fig.add_trace(go.Bar(
        x=cell_type_names,
        y=data['posterior_probs'] * 100,
        marker_color=final_colors,
        showlegend=False,
        hovertemplate='<b>%{x}</b><br>posterior: %{y:.1f}%<extra></extra>'
    ), row=4, col=1)

    # Target subplot size based on provided screenshot dimensions (426x369 px)
    target_subplot_w = 426
    target_subplot_h = 369

    # Compute overall figure size to approximate per-subplot dimensions
    # Note: Plotly spacing is fractional, so this is an approximation.
    fig_width = target_subplot_w * 2 + 160  # margins/padding
    fig_height = target_subplot_h * 4 + 240  # margins/padding

    # Update layout using computed figure size
    alpha_note = "custom" if data.get('alpha_source') == 'override' else "default"
    fig.update_layout(
        height=fig_height,
        width=fig_width,
        title_text=(
            f"<span style='font-size:18px'><b>Cell {data['label']}: Classification Trace</b></span><br>"
            f"<span style='font-size:12px'>Predicted: {data['predicted_class']} ({data['predicted_prob'] * 100:.1f}%) | "
            f"Mode: {data['prior_mode']} | Alpha: {alpha_note}</span>"
        ),
        title_x=0.5,
        title_y=0.98,
        template='plotly_white',
        font=dict(family="Arial, sans-serif", size=11),
        # Increase top margin to add padding between title and top row
        margin=dict(l=80, r=40, t=130, b=60)
    )

    # Update y-axes labels
    fig.update_yaxes(title_text="ini_alpha", row=1, col=1)
    fig.update_yaxes(title_text="ini_alpha + zeta", row=1, col=2)
    fig.update_yaxes(title_text="Log Prior", row=2, col=1)
    fig.update_yaxes(title_text="Prior (%)", row=2, col=2)
    fig.update_yaxes(title_text="Log-Likelihood", row=3, col=1)
    fig.update_yaxes(title_text="Log Posterior", row=3, col=2)
    fig.update_yaxes(title_text="Posterior (%)", row=4, col=1)

    # Update x-axes
    for row in [1, 2, 3, 4]:
        for col in [1, 2]:
            fig.update_xaxes(tickangle=-45, row=row, col=col)

    fig.show()


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
        logger.error(f"Convergence check failed: {str(e)}")
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
