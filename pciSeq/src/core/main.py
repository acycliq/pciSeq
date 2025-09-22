"""
Core Algorithm Implementation Module for pciSeq

This module implements the main Variational Bayes algorithm for spatial transcriptomics
analysis, primarily through the VarBayes class. The algorithm iteratively:
1. Assigns spots to cells
2. Determines cell types
3. Estimates gene expression patterns
4. Updates model parameters

Key Components:
-------------
VarBayes:
    Main class implementing the iterative algorithm with methods for:
    - Gene count updates
    - Cell type assignment
    - Spot-to-cell assignment
    - Parameter estimation (eta, gamma, covariance)
    - Model convergence checking

Algorithm Steps:
--------------
1. Initialization:
   - Set prior probabilities
   - Initialize cell assignments
   - Set gene efficiency parameters

2. Iterative Updates:
   - Update expected gene counts
   - Calculate gamma expectations
   - Update gaussian parameters
   - Assign cells to types
   - Assign spots to cells
   - Update gene efficiency
   - Update Dirichlet parameters
   - Update single-cell reference

3. Convergence:
   - Check for convergence after each iteration
   - Return results when converged or max iterations reached

Notes:
-----
- Uses Redis for optional diagnostic monitoring
- Implements equations from the pciSeq paper
- Handles missing single-cell reference data
- Supports parallel processing via numpy operations

Dependencies:
-----------
- numpy: For numerical computations
- pandas: For data management
- scipy: For statistical operations
- numpy_groupies: For group operations
- dask: For delayed computations
"""
import logging
from typing import Dict, List, Optional, Tuple, Union, Any

# Third-party imports
import sys
import numpy as np
import numpy_groupies as npg
import pandas as pd
from dask.delayed import delayed
from scipy.special import softmax
import opt_einsum as oe

# Local imports
from .datatypes.cells import Cells
from .datatypes.genes import Genes
from .datatypes.spots import Spots
from .datatypes.singleCell import SingleCell
from .datatypes.cellClass import CellClass
from .summary import collect_data
# from .analysis import CellExplorer
from .utils import ops_utils as utils
from .utils import visualisation
from ...src.diagnostics.controller.diagnostic_controller import DiagnosticController
import joblib

# Configure logging
main_logger = logging.getLogger(__name__)


class VarBayes:
    """
    Implements Variational Bayes algorithm for spatial transcriptomics analysis.

    This class performs cell type assignment and spot-to-cell mapping using a
    probabilistic model with variational inference.

    Args:
        cells_df: DataFrame containing cell information
        spots_df: DataFrame containing spot information
        scRNAseq: Single-cell RNA sequencing reference data
        config: Configuration dictionary containing algorithm parameters
    """

    def __init__(self,
                 cells_df: pd.DataFrame,
                 spots_df: pd.DataFrame,
                 scRNAseq: pd.DataFrame,
                 config: Dict[str, Any]) -> None:
        """Initialize components and setup."""
        # Explicitly declare important instance attributes
        self.diagnostic_controller: Optional[DiagnosticController] = None  # For real-time diagnostics
        self.config = None
        self.iter_num = None
        self.iter_delta = []
        self.has_converged = False

        # Initialize components
        self._validate_config(config)
        self.config = config
        self._setup_diagnostics()
        self._setup_components(cells_df, spots_df, scRNAseq)
        self._setup_dimensions()

        # Placeholder for other attributes
        self._scaled_exp = None
        # self._cell_explorer: Optional[CellExplorer] = None

    @staticmethod
    def _validate_config(config: Dict[str, Any]) -> None:
        """Check for required config parameters."""
        required = ['exclude_genes', 'max_iter', 'CellCallTolerance',
                    'rGene', 'Inefficiency', 'InsideCellBonus', 'MisreadDensity',
                    'cell_centroid_prior', 'cell_cov_prior', 'SpotReg', 'nNeighbors', 'rSpot',
                    'save_data', 'output_path', 'launch_viewer', 'launch_diagnostics',
                    'is_redis_running', 'cell_radius', 'cell_type_prior', 'is3D',
                    'mean_gene_counts_per_class', 'mean_gene_counts_per_cell']
        missing = [param for param in required if param not in config]
        if missing:
            raise ValueError(f"Missing required config parameters: {missing}")

    def _setup_diagnostics(self) -> None:
        """Initialize diagnostics controller if enabled in config."""
        self.diagnostic_controller = None
        if not self.config.get('launch_diagnostics', False):
            return

        try:
            self.diagnostic_controller = DiagnosticController()
            if not self.diagnostic_controller.launch_dashboard():
                main_logger.warning("Failed to launch diagnostics dashboard")
                self.diagnostic_controller = None
        except Exception as e:
            main_logger.warning(f"Failed to initialize diagnostics: {e}")
            self.diagnostic_controller = None

    def _setup_components(self, cells_df, spots_df, scRNAseq) -> None:
        """Set up the core data components needed for the algorithm."""
        self.cells = Cells(cells_df, self.config)
        self.spots = Spots(spots_df, self.config)
        self.genes = Genes(self.spots, self.config)
        self.single_cell = SingleCell(scRNAseq, self.genes.gene_panel, self.config)
        self.cellTypes = CellClass(self.single_cell, self.config)
        self.cells.class_names = self.single_cell.classes

    def _setup_dimensions(self) -> None:
        """Set up core dimensions."""
        self.nC = self.cells.nC  # cells
        self.nG = self.genes.nG  # genes
        self.nK = self.cellTypes.nK  # classes
        self.nS = self.spots.nS  # spots
        self.nN = self.config['nNeighbors'] + 1  # neighbors + background

    def initialise_state(self) -> None:
        """Initialises the starting state of the objects
        4-Feb-2025: Note that inefficiency (denoted by eta) follows a Gamma(rGene, rGene).
        Keep in mind also that there is also the config['Inefficiency'] parameter that has
        been applied directly to the expression data from scRNAseq
        """
        self.cellTypes.ini_prior()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.cells._plane_adj = self.cells.calc_plane_adj(self.spots, self.config)
        self.genes.init_eta(self.config['rGene'], self.config['rGene'])
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)[0]
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id, self.config)
        self.spots._plane_adj = self.spots.calc_plane_adj(self.cells, self.config)
        self.cells._ini_gene_counts = np.bincount(self.spots.data.label.values, minlength=self.nC)
        self.genes._misread_density = self.genes.calc_misread_density()

    def __getstate__(self):
        """
        Get state for pickling.
        Removes diagnostics-related attributes to enable pickling and reduce file size.
        """
        attributes = self.__dict__.copy()
        if 'diagnostic_controller' in attributes:
            del attributes['diagnostic_controller']
        return attributes

    @property
    def scaled_exp(self):
        """
        Get scaled expression values.

        Returns:
            delayed: Dask delayed object containing scaled expression computation
        """
        return self._scaled_exp

    # @property
    # def cell_explorer(self) -> CellExplorer:
    #     """
    #     Get cell analyzer instance.
    #     Returns:
    #         CellExplorer: Instance configured for this VarBayes object
    #     """
    #     if self._cell_explorer is None:
    #         self._cell_explorer = CellExplorer(self)
    #     return self._cell_explorer

    # -------------------------------------------------------------------- #
    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        self.initialise_state()
        cell_df, gene_df = self.main_loop()
        return cell_df, gene_df

    # -------------------------------------------------------------------- #
    def main_loop(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Main algorithm loop with diagnostic updates."""
        """
        Executes the main Variational Bayes algorithm loop.

        Iteratively updates:
            1. Gene counts per cell
            2. Gamma parameters
            3. Gaussian parameters (if needed)
            4. Cell type assignments
            5. Spot-to-cell assignments
            6. Gene efficiency parameters
            7. Dirichlet parameters (if needed)
            8. Expression means (if needed)

        The loop continues until either:
            - Convergence is reached (change in probabilities below tolerance)
            - Maximum iterations are reached

        Args:
            None

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]:
                - Cell dataframe with final assignments and probabilities
                - Gene dataframe with expression statistics

        Note:
            Progress is published to Redis if enabled.
        """
        p0 = None
        cell_df = None
        gene_df = None
        max_iter = self.config['max_iter']

        try:
            for i in range(max_iter):
                self.iter_num = i

                # 1. For each cell, calc the expected gene counts
                self.geneCount_upd()

                # 2. calc expected gamma
                self.gamma_upd()

                main_logger.info("gaussian_upd step has been removed in this version of the software")
                # 3 update correlation matrix and variance of the gaussian distribution
                # if self.single_cell.isMissing or (self.config['InsideCellBonus'] is False) or (self.config['is3D']):
                #     self.gaussian_upd()

                # 4. assign cells to cell types
                self.cell_to_cellType()

                # 5. assign spots to cells
                self.spots_to_cell()

                # 6. update gene efficiency
                self.eta_upd()

                # 7. update the dirichlet distribution
                if self.single_cell.isMissing or (self.config['cell_type_prior'] == 'weighted'):
                    self.dalpha_upd()

                # 8. Update single cell data
                if self.single_cell.isMissing:
                    self.mu_upd()

                self.has_converged, delta = utils.has_converged(
                    self.spots, p0, self.config['CellCallTolerance']
                )
                main_logger.info('Iteration %d, mean prob change %f' % (i, delta))

                # Update diagnostics using controller
                self.diagnostics_upd()

                # keep track of the deltas
                self.iter_delta.append(delta)

                # replace p0 with the latest probabilities
                p0 = self.spots.parent_cell_prob

                if self.has_converged:
                    # self.cell_analysis(35975)
                    cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.config['is3D'])
                    break

                if i == max_iter - 1:
                    main_logger.info('Loop exhausted. Exiting with convergence status: %s' % self.has_converged)
                    cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.config['is3D'])
                    break
        finally:
            # Ensure diagnostics are properly shut down
            if self.diagnostic_controller is not None:
                try:
                    self.diagnostic_controller.shutdown()
                except Exception as e:
                    main_logger.warning(f"Failed to shutdown diagnostics: {e}")

        return cell_df, gene_df

    # -------------------------------------------------------------------- #
    def geneCount_upd(self) -> None:
        """
        Updates the gene count matrix for each cell.

        Produces a matrix numCells-by-numGenes where element at position (c,g) keeps
        the expected counts of gene g in cell c. The first row corresponds to the
        background counts (spots not assigned to any cell).

        Note:
            The sum of background spots and cell gene counts should equal
            the total number of spots.
        """
        # make an array nS-by-nN and fill it with the spots id
        gene_ids = np.tile(self.spots.gene_id, (self.nN, 1)).T

        # flatten it
        gene_ids = gene_ids.ravel()

        # make corresponding arrays for cell_id and probs
        cell_ids = self.spots.parent_cell_id.ravel()
        probs = self.spots.parent_cell_prob.ravel()

        # make the array to be used as index in the group-by operation
        group_idx = np.vstack((cell_ids, gene_ids))

        # For each cell aggregate the number of spots from the same gene.
        # It will produce an array of size nC-by-nG where the entry at (c,g)
        # is the gene counts of gene g within cell c
        N_cg = npg.aggregate(group_idx, probs, size=(self.nC, self.nG))

        # assert N_cg.sum() == self.spots.data.shape[0], \
        #     "The sum of the background spots and the cell gene counts should be equal to the total number of spots"

        # make output. This part needs to be rewritten
        out = np.zeros([self.nC, self.nG], dtype=np.float32)
        out[1:, :] = N_cg[1:, :]

        # cell at position zero is the background
        self.cells.background_counts = N_cg[0, :]
        # Actual cells are on non-zero positions
        self.cells.geneCount = out

    # -------------------------------------------------------------------- #
    def gamma_upd(self) -> None:
        """
        Updates gamma parameters for the negative binomial distribution.

        Implements equation (3) of the Qian paper. Calculates the expected gamma
        values using scaled expression and spot regularization parameters.

        Updates:
            - self._scaled_exp: Delayed computation of scaled expression
            - self.spots._log_gamma_bar: Log of expected gamma values
            - self.spots._gamma_bar: Expected gamma values
        """
        cells = self.cells
        cfg = self.config

        self._scaled_exp = delayed(utils.scaled_exp(cells.ini_cell_props['area_factor'],
                                                    self.single_cell.mean_expression_adj.values))

        beta = self.scaled_exp.compute() * self.genes.eta_bar[:, None] + cfg['rSpot']

        # adjust by plane depth
        beta = np.einsum('cg,cgk->cgk', self.cells.plane_adj.values, beta)

        rho = cfg['rSpot'] + cells.geneCount

        self.spots._log_gamma_bar = delayed(self.spots.logGammaExpectation(rho, beta))
        self.spots._gamma_bar = delayed(self.spots.gammaExpectation(rho, beta))

    # -------------------------------------------------------------------- #
    def cell_to_cellType(self) -> None:
        """
        Updates cell type assignment probabilities.

        Implements equation (2) of the Qian paper. Returns an array of size
        numCells-by-numCellTypes where element in position [i,j] keeps the
        probability that cell i has cell type j.

        The computation combines:
            1. Negative binomial log-likelihood for gene expression
            2. Cell type priors
            3. Softmax normalization for final probabilities
        """

        # Get the full log-likelihood matrix using shared computation
        contr = utils.compute_gene_loglikelihood_matrix(self)

        # populate the genes' contributions to the negative loglik. Property 'nb_contr' is only useful
        # for debugging, safe to remove in the future
        self.cells.nb_contr = contr
        contr = np.sum(contr, axis=1)
        wCellClass = contr + self.cellTypes.log_prior
        pCellClass = softmax(wCellClass, axis=1)

        self.cells.classProb = pCellClass

    # -------------------------------------------------------------------- #
    def spots_to_cell(self) -> None:
        """
        Updates spot-to-cell assignment probabilities.

        Implements equation (4) of the Qian paper. For each spot, calculates the
        probability of it belonging to each nearby cell or being a misread.

        The computation includes:
            1. Expected gene expression for each cell type
            2. Gamma parameter contributions
            3. Gene efficiency factors
            4. Spatial distance likelihood
            5. Inside-cell bonus for spots within cell boundaries
            6. Misread probability for background noise

        Note:
            Updates spot-cell assignments and triggers gene count update.
        """
        nN = self.nN
        nS = self.spots.data.gene_name.shape[0]

        wSpotCell = np.zeros([nS, nN], dtype=np.float64)
        gn = self.spots.data.gene_name.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values
        logeta_bar = self.genes.logeta_bar[self.spots.gene_id]

        # misread = self.spot_misread_density()
        misread = self.spots.misread_density(self.genes)

        # pre-populate last column
        wSpotCell[:, -1] = np.log(misread)
        mvn_loglik_arr = np.zeros(wSpotCell.shape)
        attention = np.zeros(wSpotCell.shape)
        expr_fluctuations = np.zeros(wSpotCell.shape)

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # adjust the single cell data based on the plane of the cell
            spot_plane_adj = self.spots.plane_adj[:, n]  # shape: (nS,)
            expected_counts_adj = expected_counts + np.log(spot_plane_adj[:, None])  # broadcast to (nS, nK)

            # multiply and sum over cells. In practice this means that when high expected counts
            # are aligned with high cell class probs this term will be high
            term_1 = np.einsum('ij, ij -> i', expected_counts_adj, cp)

            log_gamma_bar = self.spots.log_gamma_bar.compute()
            log_gamma_bar = log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)

            # wSpotCell[:, n] = term_1 + term_2 + logeta_bar + loglik[:, n]
            mvn_loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells, self.config['is3D'])
            wSpotCell[:, n] = term_1 + term_2 + mvn_loglik
            mvn_loglik_arr[:, n] = mvn_loglik
            attention[:, n] = term_1
            expr_fluctuations[:, n] = term_2

        # apply inside cell bonus
        bonus_mask = self.spots.bonus_mask * self.config['InsideCellBonus']
        wSpotCell += bonus_mask

        # update the prob a spot belongs to a neighboring cell
        self.spots.parent_cell_prob = softmax(wSpotCell, axis=1)
        self.spots.mvn_loglik_arr = mvn_loglik_arr
        self.spots.attention = attention
        self.spots.expr_fluctuations = expr_fluctuations

        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()

    # -------------------------------------------------------------------- #
    def spots_to_cell_par(self) -> None:
        """
        Updates spot-to-cell assignment probabilities.

        Implements equation (4) of the Qian paper. For each spot, calculates the
        probability of it belonging to each nearby cell or being a misread.

        Parallelized (multithreading) version of 'spots_to_cell'
        """
        nN = self.nN
        nS = self.nS

        wSpotCell = np.zeros([nS, nN], dtype=np.float64)
        gn = self.spots.data.gene_name.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values

        # Pre-populate misread column
        misread = self.spots.misread_density(self.genes)
        wSpotCell[:, -1] = np.log(misread)

        log_gamma_bar = self.spots.log_gamma_bar.compute()

        def process_neighbor(n):
            sn = self.spots.parent_cell_id[:, n]
            cp = self.cells.classProb[sn]

            term_1 = oe.contract('ij, ij -> i', expected_counts, cp, optimize='optimal')

            current_log_gamma = log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]
            term_2 = oe.contract('ij, ij -> i', cp, current_log_gamma, optimize='optimal')

            mvn_loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells, self.config['is3D'])
            return n, term_1 + term_2 + mvn_loglik

        # Parallel processing
        results = joblib.Parallel(n_jobs=-1, backend='threading')(
            joblib.delayed(process_neighbor)(n) for n in range(nN - 1)
        )

        # Fill results back into wSpotCell
        for n, result in results:
            wSpotCell[:, n] = result

        # Apply inside cell bonus
        bonus_mask = self.spots.bonus_mask * self.config['InsideCellBonus']
        wSpotCell += bonus_mask

        # Update probabilities
        self.spots.parent_cell_prob = softmax(wSpotCell, axis=1)

        # Update gene counts
        self.geneCount_upd()

    # -------------------------------------------------------------------- #
    def eta_upd(self) -> None:
        """
        Updates gene efficiency parameters (eta).

        Implements equation (5) of the Qian paper. Calculates the expected eta values
        by combining:
            1. Total gene counts across cells
            2. Cell type probabilities
            3. Mean expression values
            4. Area factors and gamma values

        Note:
            The zero-expressing cell class is excluded from the computation.
        """
        # grand_total = self.cells.background_counts.sum() + self.cells.total_counts.sum()
        # assert round(grand_total) == self.spots.data.shape[0], \
        #     'The sum of the background spots and the total gene counts should be equal to the number of spots'

        classProb = self.cells.classProb
        mu = self.single_cell.mean_expression_adj + self.config['SpotReg']
        area_factor = self.cells.ini_cell_props['area_factor']
        gamma_bar = self.spots.gamma_bar.compute()
        plane_adj = self.cells.plane_adj.values

        zero_prob = classProb[:, -1]  # probability a cell being a zero expressing cell
        zero_class_counts = self.spots.zero_class_counts(self.spots.gene_id, zero_prob)
        # zero_class_counts = oe.contract('c, cg -> g', classProb[:, -1], self.cells.geneCount, optimize='optimal')

        # Calcs the sum in the Gamma distribution (equation 5). The zero class
        # is excluded from the sum, hence the arrays in the einsum below stop at :-1
        # Note. We should exclude the "cell" that is meant to keep the
        # misreads, ie exclude the background, hence the relevant indexing below
        # starts at 1
        class_total_counts = oe.contract('ck, gk, cg, c, cgk -> g',
                                         classProb[:, :-1],
                                         mu.values[:, :-1],
                                         plane_adj,
                                         area_factor,
                                         gamma_bar[:, :, :-1], optimize='optimal')
        # background_counts = self.cells.background_counts
        background_counts = np.bincount(self.spots.gene_id, self.spots.parent_cell_prob[:, -1], minlength=self.nG)

        # observed (ie actual) gene reads per gene
        observed = self.config['rGene'] + self.spots.counts_per_gene - background_counts - zero_class_counts

        # expected (ie predicted) gene reads per gene
        expected = self.config['rGene'] + class_total_counts

        # Finally, update gene_gamma. It will basically divide observed by expected
        # and gene inefficiency will eventually express how well a gene is detected.
        self.genes.calc_eta(observed, expected)

    # -------------------------------------------------------------------- #
    def gaussian_upd(self) -> None:
        """
        Updates Gaussian distribution (centroids and covariance matrices) for cells
        """
        self.centroid_upd()
        self.cov_upd()

    # -------------------------------------------------------------------- #
    def centroid_upd(self) -> None:
        """
        Updates the centroid (mean) for each cell based on the posterior distribution.

        The posterior centroid is calculated as a weighted average of:
        - The prior centroid (mu_0), weighted by the prior pseudo-sample size (k_0)
        - The empirical (sample) mean (x_bar), weighted by the observed sample size (n)
        """

        # Get the prior weight (pseudo-sample size for the centroid)
        k_0 = self.config['cell_centroid_prior']['default']

        # Prior centroid (mu_0)
        prior_centroid = self.cells.ini_centroids()

        # 1. Calculate the empirical (sample) mean
        sample_mean = utils.empirical_mean(spots=self.spots, cells=self.cells)

        # 2. Get the observed sample size (gene counts per cell)
        sample_size = self.cells.total_counts

        # 3. Calculate the posterior centroid as a weighted average
        numerator = k_0 * prior_centroid + sample_size[:, None] * sample_mean
        denominator = k_0 + sample_size[:, None]
        posterior_centroid = numerator / denominator

        # 4. Handle the background (index 0)
        # Reset the background centroid to a default large value (e.g., max int)
        posterior_centroid.iloc[0, :] = -sys.maxsize

        # Update the centroids
        self.cells.centroid = posterior_centroid

    # -------------------------------------------------------------------- #
    def cov_upd(self) -> None:
        """
        Updates the covariance matrix for each cell based on the posterior distribution.

        The posterior covariance is computed using:
        - The scatter matrix (data-driven scale matrix)
        - The prior covariance (weighted by prior hyperparameters)
        - An adjustment term based on the difference between sample means and prior means
        """

        # Hyperparameters
        k_0 = self.config['cell_centroid_prior']['default']  # Prior for centroids
        nu_0 = self.config['cell_cov_prior']['default']  # Prior degrees of freedom for covariance

        # 1. Calculate the scatter matrix (data-driven scale matrix)
        scatter_matrix = self.cells.scatter_matrix(self.spots)

        # 2. Calculate the prior scale matrix
        prior_cov = self.cells.ini_cov()
        prior_scale_matrix = nu_0 * prior_cov

        # 3. Calculate the adjustment term
        # Difference between current centroids and prior centroids (x_bar - mu_0)
        mean_diff = self.cells.centroid - self.cells.ini_centroids()
        mean_outer_product = oe.contract('rk, rn -> rkn', mean_diff, mean_diff,
                                         optimize='optimal')  # (x_bar - mu_0)(x_bar - mu_0)^T

        # Multiplier for the adjustment term
        multiplier = (k_0 * self.cells.total_counts) / (k_0 + self.cells.total_counts)

        # Avoid warnings by setting background (index 0) adjustment to zero
        mean_outer_product[0] = np.zeros_like(mean_outer_product[0])
        adjustment_term = multiplier[:, None, None] * mean_outer_product

        # 4. Calculate the updated scale matrix
        scale_matrix_upd = scatter_matrix + prior_scale_matrix + adjustment_term

        # 5. Calculate the updated degrees of freedom
        nu_upd = nu_0 + self.cells.total_counts + 1

        # 6. Compute the expected covariance
        covariance_upd = utils.expected_covariance(scale_matrix_upd, nu_upd)

        # Handle background (index 0) by mapping it to the prior covariance
        covariance_upd[0] = prior_cov[0]

        # 7. Update the cell covariance attribute. Get the eigenvals/eigenvectors too
        self.cells.cov = covariance_upd
        self.cells.eig_vals, self.cells.eig_vecs = np.linalg.eigh(covariance_upd)

    # -------------------------------------------------------------------- #
    def mu_upd(self) -> None:
        """
        Updates mean expression values when single-cell reference is missing.

        Estimates mean expression for each gene and cell type using:
            1. Current cell type assignments
            2. Observed gene counts
            3. Cell area factors
            4. Current gamma and eta values

        Updates:
            - single_cell._mean_expression: Updated mean expression values
            - single_cell._log_mean_expression: Log of mean expression values
        """
        classProb = self.cells.classProb[1:, :-1].copy()
        geneCount = self.cells.geneCount[1:, :].copy()
        gamma_bar = self.spots.gamma_bar.compute()[1:, :, :-1]
        area_factor = self.cells.ini_cell_props['area_factor'][1:]

        numer = oe.contract('ck, cg -> gk', classProb, geneCount, optimize='optimal')
        denom = oe.contract('ck, c, cgk, g -> gk', classProb, area_factor, gamma_bar, self.genes.eta_bar,
                            optimize='optimal')

        me, lme = self.single_cell._gene_expressions(numer, denom)
        self.single_cell._mean_expression = me
        self.single_cell._log_mean_expression = lme

    # -------------------------------------------------------------------- #
    def dalpha_upd(self) -> None:
        """
        Updates cell type prior distribution parameters.

        Adjusts Dirichlet parameters based on:
            1. Current cell type assignments
            2. Initial alpha values
            3. Minimum class size constraints

        Note:
            - Ensures Zero class (background) is preserved
            - Sets very small weights (1e-6) for classes below minimum size
        """
        # logger.info('Update cell type (marginal) distribution')
        zeta = self.cells.classProb.sum(axis=0)  # this the class size
        alpha = self.cellTypes.ini_alpha()
        out = zeta + alpha

        # 07-May-2023: Hiding 'min_class_size' from the config file. Should bring it back at a later version
        # mask = zeta <= self.config['min_class_size']
        min_class_size = 5
        mask = zeta <= min_class_size

        # make sure Zero class is the last one
        assert self.cellTypes.names[-1] == "Zero"
        assert len(self.cellTypes.names) == len(mask)

        # make sure the last value ie always False, overriding if necessary the
        # check a few lines above when the mask variable was set.
        # In this manner we will prevent the Zero class from being removed.
        mask[-1] = False

        # If a class size is smaller than 'min_class_size' then it will be assigned a weight of almost zero
        out[mask] = 10e-6
        self.cellTypes.alpha = out

    # -------------------------------------------------------------------- #
    # def spot_misread_density(self) -> np.array:
    #     """
    #     Calculates spot misread probabilities for each gene.
    #
    #     Combines:
    #         1. Default misread probability for all genes
    #         2. Gene-specific probabilities from configuration
    #         3. Alignment with current spot assignments
    #
    #     Returns:
    #         np.ndarray: Array of misread probabilities aligned with spots
    #     """
    #     # Get default misread probability for all genes
    #     default_val = self.config['MisreadDensity']['default']
    #     gene_names = self.genes.gene_panel
    #     misread_dict = dict(zip(gene_names, [default_val] * self.nG))
    #
    #     # Update with any gene-specific probabilities
    #     misread_dict.update(self.config['MisreadDensity'] or {})
    #     misread_dict.pop('default', None)
    #
    #     # Convert to array and align directly with spots
    #     v = np.array(list(misread_dict.values()))
    #     v = v[self.spots.gene_id]  # Align with spots
    #     return v

    # -------------------------------------------------------------------- #
    def diagnostics_upd(self) -> None:
        """Update diagnostic visualization if controller is available."""
        if self.diagnostic_controller is None:
            return

        try:
            self.diagnostic_controller.update_diagnostics(
                algorithm_model=self,
                iteration=self.iter_num,
                has_converged=self.has_converged
            )
        except Exception as e:
            main_logger.warning(f"Failed to update diagnostics: {e}")

    # -------------------------------------------------------------------- #
    # def cell_analysis(self, cell_num):
    #     """
    #     Convenience method to analyze a specific cell.
    #
    #     Parameters
    #     ----------
    #     cell_num : int
    #         The cell number to analyze
    #
    #     Returns
    #     -------
    #     Same as cell_explorer.view_cell()
    #     """
    #     return self.cell_explorer.view_cell(cell_num)

    # -------------------------------------------------------------------- #
    def heatmap_counts_per_class(self):
        """Display the interactive heatmap."""
        return visualisation.heatmap_counts_per_class(self)

    def calculate_genes_log_likelihood_contr(self, label):
        return utils.calculate_genes_log_likelihood_contr(self, label)

    # def plot_loglik_contr(self, df):
    #     return utils.plot_loglik_contr(df)

    # def visualize_fit(self, gene_counts, scaled_means):
    #     return utils.visualize_fit(gene_counts, scaled_means)

    def check_cell(self, my_label, user_class, top_n=10, show_plot=True):
        return utils.check_cell(self, my_label, user_class, top_n, show_plot)

    def check_spot(self, spot_id):
        return visualisation.check_spot(self, spot_id)

    def read_tsv(self, filepath):
        return utils.read_tsv(filepath)

    # def trellis_plot(self, label, flatfile_folder):
    #     return visualisation.trellis_plot(self, label, flatfile_folder)


