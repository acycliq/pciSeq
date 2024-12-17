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
import numpy as np
import numpy_groupies as npg
import pandas as pd
from dask.delayed import delayed
import scipy.spatial as spatial
from scipy.special import softmax

# Local imports
from .datatypes import Cells, Spots, Genes, SingleCell, CellType
from .summary import collect_data
from .analysis import CellExplorer
from .utils import ops_utils as utils
from ...src.diagnostics.controller.diagnostic_controller import DiagnosticController

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
        self._cell_explorer: Optional[CellExplorer] = None

    @staticmethod
    def _validate_config(config: Dict[str, Any]) -> None:
        """Check for required config parameters."""
        required = ['exclude_genes', 'max_iter', 'CellCallTolerance',
                    'rGene', 'Inefficiency', 'InsideCellBonus', 'MisreadDensity',
                    'cell_centroid_prior_weight', 'SpotReg', 'nNeighbors', 'rSpot',
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
        self.genes = Genes(self.spots)
        self.single_cell = SingleCell(scRNAseq, self.genes.gene_panel, self.config)
        self.cellTypes = CellType(self.single_cell, self.config)
        self.cells.class_names = self.single_cell.classes

    def _setup_dimensions(self) -> None:
        """Set up core dimensions."""
        self.nC = self.cells.nC  # cells
        self.nG = self.genes.nG  # genes
        self.nK = self.cellTypes.nK  # classes
        self.nS = self.spots.nS  # spots
        self.nN = self.config['nNeighbors'] + 1  # neighbors + background

    def initialise_state(self) -> None:
        self.cellTypes.ini_prior()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.genes.init_eta(1, 1 / self.config['Inefficiency'])
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)[0]
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id, self.config)

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

    @property
    def cell_explorer(self) -> CellExplorer:
        """
        Get cell analyzer instance.
        Returns:
            CellExplorer: Instance configured for this VarBayes object
        """
        if self._cell_explorer is None:
            self._cell_explorer = CellExplorer(self)
        return self._cell_explorer

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

        self.initialise_state()
        try:
            for i in range(max_iter):
                self.iter_num = i

                # 1. For each cell, calc the expected gene counts
                self.geneCount_upd()

                # 2. calc expected gamma
                self.gamma_upd()

                # 3 update correlation matrix and variance of the gaussian distribution
                if self.single_cell.isMissing or (self.config['InsideCellBonus'] is False) or (self.config['is3D']):
                    self.gaussian_upd()

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
                    # self.cell_analysis(2259)
                    cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell,
                                                    self.config['is3D'])
                    break

                if i == max_iter - 1:
                    main_logger.info('Loop exhausted. Exiting with convergence status: %s' % self.has_converged)

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
                                                    self.single_cell.mean_expression.values,
                                                    self.genes.eta_bar))

        beta = self.scaled_exp.compute() + cfg['rSpot']
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

        ScaledExp = self.scaled_exp.compute()
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negative_binomial_loglikelihood(cgc, self.config['rSpot'], pNegBin)
        contr = np.sum(contr, axis=1)
        wCellClass = contr + self.cellTypes.log_prior
        pCellClass = softmax(wCellClass, axis=1)
        del contr

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

        misread = self.spot_misread_density()

        # pre-populate last column
        wSpotCell[:, -1] = np.log(misread)

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            term_1 = np.einsum('ij, ij -> i', expected_counts, cp)

            log_gamma_bar = self.spots.log_gamma_bar.compute()
            log_gamma_bar = log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)

            # wSpotCell[:, n] = term_1 + term_2 + logeta_bar + loglik[:, n]
            mvn_loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells, self.config['is3D'])
            wSpotCell[:, n] = term_1 + term_2 + logeta_bar + mvn_loglik
            del term_1

        # apply inside cell bonus
        # NOTE. This is not applied 100% correctly. For example the fourth spot in the demo data. Id2, (x, y) = (0, 4484)
        # The spots is within the boundaries of cell with label 5 but this is not its closest cell. The closest cell is
        # cell label = 4. Cell label=5 is the second closest cell and label=4 the first closest. Therefore, the bonus should be
        # applied when we handle the column for the seconds closest near-by cell. The implementation below implies that if
        # a spot is inside the cell boundaries then that cell is the closest one.
        mask = np.greater(self.spots.data.label.values, 0, where=~np.isnan(self.spots.data.label.values))
        wSpotCell[mask, 0] = wSpotCell[mask, 0] + self.config['InsideCellBonus']

        # update the prob a spot belongs to a neighboring cell
        self.spots.parent_cell_prob = softmax(wSpotCell, axis=1)

        # Since the spot-to-cell assignments changed you need to update the gene counts now
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
        grand_total = self.cells.background_counts.sum() + self.cells.total_counts.sum()
        assert round(grand_total) == self.spots.data.shape[0], \
            'The sum of the background spots and the total gene counts should be equal to the number of spots'

        classProb = self.cells.classProb
        mu = self.single_cell.mean_expression
        area_factor = self.cells.ini_cell_props['area_factor']
        gamma_bar = self.spots.gamma_bar.compute()

        zero_prob = classProb[:, -1]  # probability a cell being a zero expressing cell
        zero_class_counts = self.spots.zero_class_counts(self.spots.gene_id, zero_prob)

        # Calcs the sum in the Gamma distribution (equation 5). The zero class
        # is excluded from the sum, hence the arrays in the einsum below stop at :-1
        # Note. I also think we should exclude the "cell" that is meant to keep the
        # misreads, ie exclude the background
        class_total_counts = np.einsum('ck, gk, c, cgk -> g',
                                       classProb[:, :-1],
                                       mu.values[:, :-1],
                                       area_factor,
                                       gamma_bar[:, :, :-1])
        background_counts = self.cells.background_counts
        alpha = self.config['rGene'] + self.spots.counts_per_gene - background_counts - zero_class_counts
        beta = self.config['rGene'] / self.config['Inefficiency'] + class_total_counts

        # Finally, update gene_gamma
        self.genes.calc_eta(alpha, beta)

    # -------------------------------------------------------------------- #
    def gaussian_upd(self) -> None:
        """
        Updates Gaussian distribution (centroids and covariance matrices) for cells
        """
        self.centroid_upd()
        self.cov_upd()

    # -------------------------------------------------------------------- #
    def centroid_upd(self):
        """
        Updates cell centroids using Bayesian posterior estimation.

        Combines prior centroids with empirical means using configurable weights:
        mu_post = (alpha * mu_0 + x_bar) / (alpha + 1)

        where:
        - mu_0 is the prior centroid
        - x_bar is the empirical mean
        - alpha is the weight controlling prior influence

        Technical Details:
        -----------------
        Computes the posterior cluster mean given prior and sample statistics.

        The posterior cluster mean is derived under the assumption of a Normal-Inverse Wishart (NIW) prior, where
        the cluster's mean mu and covariance Sigma are jointly modeled.

        When given sample data consisting of `n` observations with a sample mean x_bar, the posterior mean mu_n is
        computed as:

            mu_post = (k_0 * mu_0 + n * x_bar) / (k_0 + n)

        The parameter k_0 controls the weight of the prior mean relative to the data when calculating the posterior
        mean.
         - When k_0 = 0 then the posterior mean is the sample mean
         - A large value of k_0 implies stronger confidence in the prior, making the posterior mean closer to mu_0

        Alternative Formula:
        k_0 as a Function of alpha, ie: k_0 = alpha * n
        ------------------------------------------------
        k_0 can be set proportional to the sample size n as k_0 = alpha * n, where alpha >= 0 and controls the
        relative influence of the prior:

            mu_post = (alpha * mu_0 + x_bar) / (alpha + 1)

        Setting alpha:
        ----------------
        - alpha = 0: Fully data-driven, no prior information.
        - alpha = 1: Balance the prior and the data equally, mu_post = (mu_0 + x_bar) / 2.
        - alpha >> 1: Strong prior influence, posterior mean gets closer to the prior mean.
        - alpha -> infinity: Prior dominates, posterior mean appx equal to the prior mean.
        """

        # Get default value for the weight
        default_val = self.config['cell_centroid_prior_weight']['default']

        cell_labels = np.arange(self.nC)
        alpha_dict = {label: default_val for label in cell_labels}

        # Update with any cell specific weights (Review this, can be done in a simpler way!)
        if self.config['label_map']:
            _label = []
            for key in self.config['cell_centroid_prior_weight'].keys():
                if key == 'default':
                    _label.append('default')
                else:
                    try:
                        _label.append(self.config['label_map'][key])
                    except KeyError as e:
                        main_logger.warning(f"Could not find cell with label: {key}, "
                                            f"cell_centroid_prior_weight is not applied")
            _dict = dict(zip(_label, self.config['cell_centroid_prior_weight'].values()))
        else:
            _dict = self.config['cell_centroid_prior_weight']

        alpha_dict.update(_dict or {})
        alpha_dict.pop('default', None)

        mu_0 = self.cells.ini_centroids()

        # 1. calc the empirical (sample) mean
        x_bar = utils.empirical_mean(spots=self.spots, cells=self.cells)

        # 2.  Weighted average of the prior mean and the empirical mean
        alpha = np.fromiter(alpha_dict.values(), dtype=np.float32)
        alpha = alpha[:, None]
        a = alpha * mu_0 + x_bar
        b = alpha + 1
        mu_post = a / b
        self.cells.centroid = mu_post

    # -------------------------------------------------------------------- #
    def cov_upd(self) -> None:
        """
        TBA
        nu_0 is the degrees of freedom. This reflects how much trust you place in your prior belief
        about the covariance structure before seeing any data.
        A larger nu_0 means you trust the prior covariance structure more because the prior is
        "informed" by more data or is based on a more confident belief.
        For example, in the context of a covariance matrix, \nu_0 could represent how much "prior knowledge"
        or prior data you have regarding the covariance before observing new data.

        Impact of Degrees of Freedom:
        Small \nu_0
        (few prior data points or weak prior knowledge): The posterior estimate will rely more heavily on
        the data, and the posterior distribution will be closer to the sample covariance matrix

        Large \nu_0
        (strong prior knowledge): The prior covariance matrix will have a larger influence on the posterior
        distribution, and the posterior covariance matrix will be more influenced by the prior structure.

        In summary, degrees of freedom essentially determine how much weight is given to the prior information
        versus the data when estimating parameters like covariance.
        """
        spots = self.spots
        n = self.cells.geneCount.sum(axis=1)  # sample size (cell gene counts)
        k_0 = 20  # maybe set this equal to degrees of freedom, nu_0
        d = 3 if self.config['is3D'] else 2  # dimensionality of the data points

        # 1. first set the prior scale matrix
        cov_0 = self.cells.ini_cov()
        nu_0 = self.cells.nu_0
        psi_0 = cov_0 * nu_0

        # 2. Get now the scatter matrix. This is basically the sample covariance matrix
        # scaled be sample size (cell gene counts)
        S = self.cells.scatter_matrix(spots)

        a = S + cov_0 * nu_0
        b = n + nu_0 - d - 1

        # divide a by b (same as a/b[:, :, None])
        cov = np.einsum('crk, c -> crk', a, 1 / b)
        self.cells.cov = cov.astype(np.float32)

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

        numer = np.einsum('ck, cg -> gk', classProb, geneCount)
        denom = np.einsum('ck, c, cgk, g -> gk', classProb, area_factor, gamma_bar, self.genes.eta_bar)

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
    def spot_misread_density(self) -> np.array:
        """
        Calculates spot misread probabilities for each gene.

        Combines:
            1. Default misread probability for all genes
            2. Gene-specific probabilities from configuration
            3. Alignment with current spot assignments

        Returns:
            np.ndarray: Array of misread probabilities aligned with spots
        """
        # Get default misread probability for all genes
        default_val = self.config['MisreadDensity']['default']
        gene_names = self.genes.gene_panel
        misread_dict = dict(zip(gene_names, [default_val] * self.nG))

        # Update with any gene-specific probabilities
        misread_dict.update(self.config['MisreadDensity'] or {})
        misread_dict.pop('default', None)

        # Convert to array and align directly with spots
        v = np.array(list(misread_dict.values()))
        v = v[self.spots.gene_id]  # Align with spots
        return v

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
    def cell_analysis(self, cell_num):
        """
        Convenience method to analyze a specific cell.

        Parameters
        ----------
        cell_num : int
            The cell number to analyze

        Returns
        -------
        Same as cell_explorer.view_cell()
        """
        return self.cell_explorer.view_cell(cell_num)
