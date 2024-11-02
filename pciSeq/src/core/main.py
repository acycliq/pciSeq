# Standard library imports
import json
import logging
import os
import shutil
import time
from typing import Dict, List, Optional, Tuple, Union

# Third-party imports
import numpy as np
import numpy_groupies as npg
import pandas as pd
import scipy.spatial as spatial
from dask.delayed import delayed
from scipy.special import softmax
from dataclasses import dataclass
from typing import Dict, Any

# Local imports
from pciSeq.src.core.datatypes import Cells, Spots, Genes, SingleCell, CellType
from pciSeq.src.core.summary import collect_data
from pciSeq.src.core import utils
from pciSeq.src.diagnostics.redis_publisher import RedisPublisher
from pciSeq.src.diagnostics.utils import RedisDB

# Configure logging
main_logger = logging.getLogger(__name__)


@dataclass
class RedisConfig:
    """Data class for Redis configuration."""
    enabled: bool
    host: str = 'localhost'
    port: int = 6379
    db: int = 0
    flush: bool = True


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
        self._validate_config(config)
        self.config = config
        self._setup_redis()
        self._setup_components(cells_df, spots_df, scRNAseq)
        self._setup_dimensions()

        # Initialize algorithm state tracking
        self.iter_num = None
        self.iter_delta = []
        self.has_converged = False

        # Will hold computed expression values
        self._scaled_exp = None

    def _validate_config(self, config: Dict[str, Any]) -> None:
        """Check for required config parameters."""
        required = ['exclude_genes', 'max_iter', 'CellCallTolerance',
                    'rGene', 'Inefficiency', 'InsideCellBonus', 'MisreadDensity',
                    'SpotReg', 'nNeighbors', 'rSpot', 'save_data', 'output_path',
                    'launch_viewer', 'launch_diagnostics', 'is_redis_running',
                    'cell_radius', 'cell_type_prior', 'is3D', 'mean_gene_counts_per_class',
                    'mean_gene_counts_per_cell']
        missing = [param for param in required if param not in config]
        if missing:
            raise ValueError(f"Missing required config parameters: {missing}")

    def _setup_redis(self) -> None:
        """Configure Redis if enabled, otherwise set to None."""
        if self.config.get('is_redis_running', False):
            self.redis_db = RedisDB(flush=True)
            self.redis_publisher = RedisPublisher(self.redis_db)
            main_logger.info("Redis connection established")
        else:
            self.redis_db = None
            self.redis_publisher = None

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
        Removes Redis-related attributes to enable pickling and reduce file size.

        Returns:
            Dict: Object state dictionary without Redis attributes
        """
        attributes = self.__dict__.copy()
        del attributes['redis_db']
        del attributes['redis_publisher']
        # del attributes['_scaled_exp']
        return attributes

    @property
    def scaled_exp(self):
        """
        Get scaled expression values.

        Returns:
            delayed: Dask delayed object containing scaled expression computation
        """
        return self._scaled_exp

    # -------------------------------------------------------------------- #
    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        self.initialise_state()
        cell_df, gene_df = self.main_loop()
        return cell_df, gene_df

    # -------------------------------------------------------------------- #
    def main_loop(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
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

            self.has_converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            main_logger.info('Iteration %d, mean prob change %f' % (i, delta))

            # keep track of the deltas
            self.iter_delta.append(delta)

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            if self.redis_publisher:
                self.redis_publisher.publish_diagnostics(self, self.iter_num, self.has_converged)

            if self.has_converged:
                # self.counts_within_radius(20)
                # self.cell_analysis(2259)
                cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell, self.config['is3D'])
                break

            if i == max_iter - 1:
                main_logger.info('Loop exhausted. Exiting with convergence status: %s' % self.has_converged)
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
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
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
        Updates Gaussian distribution parameters for spatial modeling.

        Updates both centroids and covariance matrices for cells based on:
            1. Current spot assignments
            2. Spatial coordinates
            3. Prior parameters

        This method is called when:
            - Single cell data is missing
            - Inside cell bonus is disabled
            - 3D analysis is being performed
        """
        self.centroid_upd()
        self.cov_upd()

    # -------------------------------------------------------------------- #
    def centroid_upd(self) -> None:
        """
        Updates cell centroid positions based on spot assignments.

        Calculates weighted average positions of spots assigned to each cell.
        For cells with no assigned spots, maintains initial centroid positions.

        Updates:
            self.cells.centroid: DataFrame with updated x, y, z coordinates
        """
        spots = self.spots

        # get the total gene counts per cell
        N_c = self.cells.total_counts

        xyz_spots = spots.xyz_coords
        prob = spots.parent_cell_prob
        n = self.cells.config['nNeighbors'] + 1

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
        ini_cent = self.cells.ini_centroids()
        xyz_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y'], ini_cent['z']])))

        # # sanity check. NaNs or Infs should appear together
        # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
        # use the fitted centroids where possible otherwise use the initial ones
        xyz_bar[np.isfinite(x_bar)] = xyz_bar_fitted[np.isfinite(x_bar)]
        self.cells.centroid = pd.DataFrame(xyz_bar, columns=['x', 'y', 'z'], dtype=np.float32)
        # print(np.array(list(zip(x_bar.T, y_bar.T))))

    # -------------------------------------------------------------------- #
    def cov_upd(self) -> None:
        """
        Updates cell covariance matrices for spatial distributions.

        Combines:
            1. Empirical scatter matrix from assigned spots
            2. Prior covariance scaled by prior degrees of freedom
            3. Shrinkage to control estimation stability

        Note:
            Uses delta=0.5 for shrinkage towards prior covariance.
        """
        spots = self.spots

        # first get the scatter matrix
        S = self.cells.scatter_matrix(spots)  # sample sum of squares
        cov_0 = self.cells.ini_cov()
        nu_0 = self.cells.nu_0
        S_0 = cov_0 * nu_0  # prior sum of squares
        N_c = self.cells.total_counts
        d = 2
        denom = np.ones([self.nC, 1, 1])
        denom[:, 0, 0] = N_c + nu_0
        # Note: need to add code to handle the case N_c + nu_0 <= d + 2
        cov = (S + S_0) / denom

        # shrinkage
        delta = 0.5
        cov = delta * cov_0 + (1 - delta) * cov
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
    def counts_within_radius(self, r) -> None:
        """
        Calculates gene counts within specified radius of cell centroids.

        Args:
            r: Radius for counting spots (same units as spot coordinates)

        Returns:
            pd.DataFrame: DataFrame containing:
                - cell_label: Cell identifier
                - cell_type: Assigned cell type
                - cell_label_old: Original cell label (if relabeled)
                - Gene counts columns for each gene

        Note:
            Excludes background (label=0) from the output.
        """

        # check that the background (label=0) is at the top row
        assert self.cells.ini_cell_props['cell_label'][0] == 0

        # spots = pd.concat([self.spots.data, self.spots.data_excluded])
        gene_names, gene_id = np.unique(self.spots.data.gene_name.values,
                                        return_inverse=True)  # that needs to be encapsulated!!! Make a function in the class to set the ids
        spots = self.spots.data.assign(gene_id=gene_id)

        xy_coords = spots[['x', 'y']].values
        point_tree = spatial.cKDTree(xy_coords)
        nearby_spots = point_tree.query_ball_point(self.cells.centroid, r)

        out = np.zeros([self.cells.centroid.shape[0], len(gene_names)])

        for i, d in enumerate(nearby_spots):
            t = spots.gene_id[d]
            b = np.bincount(t)
            if len(gene_names) - len(b) < 0:
                print('oops')
            out[i, :] = np.pad(b, (0, len(gene_names) - len(b)), 'constant')

        # for each cell get the most likely cell type
        cell_type = []
        for i, d in enumerate(self.cells.classProb):
            j = self.cells.class_names[np.argmax(d)]
            cell_type.append(j)

        temp = pd.DataFrame({
            'cell_label': self.cells.ini_cell_props['cell_label'],
            'cell_type': cell_type
        })

        # if labels have been reassigned, there should be a column 'cell_label_old'.
        # Relabelling will happen if for example the original labelling that was
        # assigned from segmentation is not a continuous sequence of integers
        # Append the original labels to the temp folder so there is a mapping between
        # old and new labels
        if 'cell_label_old' in self.cells.ini_cell_props.keys():
            temp['cell_label_old'] = self.cells.ini_cell_props['cell_label_old']
            temp = temp[['cell_label', 'cell_label_old', 'cell_type']]

        # assert np.all(temp.cell_label==out.index)
        df = pd.DataFrame(out, index=self.cells.centroid.index, columns=gene_names)
        df = pd.merge(temp, df, right_index=True, left_on='cell_label', how='left')

        # ignore the first row, it is the background
        return df.iloc[1:, :]

    # -------------------------------------------------------------------- #
    def gene_loglik_contributions(self, cell_num, user_class=None) -> Dict:
        """
        Calculate and return gene log-likelihood contributions for a specified cell.

        This function analyzes how each gene contributes to the cell type classification
        by calculating log-likelihood values for each gene under different cell type hypotheses.

        Args:
            cell_num (int): The cell number to analyze. Must be between 0 and nC-1.
            user_class (str, optional): The cell class to compare against the assigned class.
                If None, uses the assigned class.

        Returns:
            dict: A dictionary containing:
                - assigned_class (str): The automatically assigned cell class
                - user_class (str): The user-specified class for comparison
                - assigned_contr (list): Log-likelihood contributions for assigned class
                - cell_num (int): The analyzed cell number
                - gene_names (list): List of gene names
                - class_names (list): List of available cell type classes
                - class_probs (dict): Probability distribution over cell types
                - contr (dict): Log-likelihood contributions for all classes

        Raises:
            ValueError: If cell_num is invalid or user_class is not recognized
            """

        if cell_num < 0 or cell_num >= self.nC:
            raise ValueError(f"Invalid cell number. Must be between 0 and {self.nC - 1}")

        assigned_class_idx = np.argmax(self.cells.classProb[cell_num])
        assigned_class = self.cellTypes.names[assigned_class_idx]

        if user_class is None:
            user_class = assigned_class

        try:
            user_class_idx = np.where(self.cellTypes.names == user_class)[0][0]
        except IndexError:
            raise ValueError(
                f"Invalid user class: {user_class}. Available classes are: {', '.join(self.cellTypes.names)}")

        # it will be nice to avoid duplicating the code. This has already been calculated.
        ScaledExp = self.scaled_exp.compute()
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)

        # Calculate contributions for all classes
        all_class_contrs = contr[cell_num, :, :]

        # Prepare the user_data dictionary with contributions for all classes
        user_data = {
            class_name: all_class_contrs[:, class_idx].tolist()
            for class_idx, class_name in enumerate(self.cellTypes.names)
        }

        class_probs = dict(zip(self.cellTypes.names.tolist(), self.cells.classProb[cell_num].tolist()))

        out = {
            'assigned_class': assigned_class,
            'user_class': user_class,
            'assigned_contr': all_class_contrs[:, assigned_class_idx].tolist(),
            'cell_num': cell_num,
            'gene_names': self.genes.gene_panel.tolist(),
            'class_names': self.cellTypes.names.tolist(),
            'class_probs': class_probs,
            'contr': user_data
        }

        # Call the plotting function
        # utils.gene_loglik_contributions_scatter(out)

        return out

    def spot_dist_and_prob(self, cell_num) -> Dict:
        """
        Calculate the relationship between spot-to-cell distances and their assignment
        probabilities for a given cell.

        This function analyzes spatial relationships between spots and a target cell by:
        1. Finding spots near the target cell
        2. Calculating distances from spots to cell centroid
        3. Computing assignment probabilities

        Args:
            cell_num (int): The cell number to analyze

        Returns:
            dict: A dictionary containing plot data:
                - x (list): Distances from spots to cell centroid
                - y (list): Assignment probabilities
                - labels (list): Gene names for each spot
                - cell_num (int): The analyzed cell number
                - title (str): Plot title
                - xlabel (str): X-axis label
                - ylabel (str): Y-axis label

        Note:
            The returned data is structured for visualization in the cell analysis dashboard.
        """

        # Get cell centroid
        centroid_zyx = self.cells.zyx_coords[cell_num]

        # Find spots near target cell
        is_spot_near_target_cell = self.spots.parent_cell_id == cell_num
        mask = np.any(is_spot_near_target_cell, axis=1)

        # Get probabilities for spots assigned to this cell
        prob = self.spots.parent_cell_prob[is_spot_near_target_cell]

        # Select relevant spots
        spots = self.spots.data[mask]

        # Find most likely parent cell for each spot
        max_idx = np.argmax(self.spots.parent_cell_prob[mask], axis=1)
        assigned_cell = np.choose(max_idx, self.spots.parent_cell_id[mask].T)

        # Calculate distances and create DataFrame
        spots = spots.assign(
            cell_x=centroid_zyx[2],
            cell_y=centroid_zyx[1],
            cell_z=centroid_zyx[0],
            prob=prob,
            assigned_cell=assigned_cell,
            dist=np.sqrt(
                (spots.x - centroid_zyx[2]) ** 2 +
                (spots.y - centroid_zyx[1]) ** 2 +
                (spots.z - centroid_zyx[0]) ** 2
            )
        )

        # Select and order columns
        spots = spots[
            ['x', 'y', 'z', 'gene_name', 'cell_x', 'cell_y', 'cell_z',
             'prob', 'assigned_cell', 'dist']
        ].reset_index()

        # Prepare plot data
        # Prepare plot data with cell-specific axis labels
        data = {
            'x': spots.dist.tolist(),
            'y': spots.prob.tolist(),
            'labels': spots.gene_name.tolist(),
            'cell_num': cell_num,
            'title': f'Cell {cell_num} - Distance vs Assignment Probability',
            'xlabel': f'Distance from cell {cell_num} centroid',
            'ylabel': f'Assignment probability to cell {cell_num}'
        }
        return data

    def cell_analysis(self, cell_num, output_dir=None) -> None:
        """
        Generates data and launches the cell analysis dashboard for a specific cell.

        This function:
        1. Generates analysis data for the specified cell
        2. Saves data to JSON files
        3. Sets up and launches a local server
        4. Opens the analysis dashboard in a web browser

        Args:
            cell_num: The cell number to analyze
            output_dir: Optional directory to save JSON files. Defaults to cell_analysis/

        Note:
            Creates a local server to serve the dashboard. Close terminal to stop server.
        """
        # Get default output directory if none specified
        if output_dir is None:
            output_dir = utils.get_out_dir(self.config['output_path'])
            output_dir = os.path.join(output_dir, 'debug', 'cell_analysis')

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Get the data2
        assigned_class_idx = np.argmax(self.cells.classProb[cell_num])
        user_class = self.cellTypes.names[assigned_class_idx]

        # Generate gene contribution data
        spot_dist = self.spot_dist_and_prob(cell_num)

        loglik_data = self.gene_loglik_contributions(cell_num, user_class)
        with open(os.path.join(output_dir, 'gene_loglik_contr.json'), 'w') as fp:
            json.dump(loglik_data, fp)
            main_logger.info(f'saved at {os.path.join(output_dir, "gene_loglik_contr.json")}')

        # Save the data files
        with open(os.path.join(output_dir, "spot_dist.json"), "w") as f:
            json.dump(spot_dist, f)
            main_logger.info(f'saved at {os.path.join(output_dir, "spot_dist.json")}')

        pciSeq_dir = utils.get_pciSeq_install_dir()
        src = os.path.join(pciSeq_dir, 'static', 'cell_analysis')

        shutil.copytree(src, output_dir, dirs_exist_ok=True)
        # viewer_utils_logger.info('viewer code (%s) copied from %s to %s' % (dim, src, dst))

        # Launch the dashboard
        import webbrowser
        import http.server
        import socketserver
        import threading
        import random

        # Start HTTP server
        os.chdir(output_dir)
        PORT = 8000 + random.randint(0, 999)

        Handler = http.server.SimpleHTTPRequestHandler

        def start_server():
            with socketserver.TCPServer(("", PORT), Handler) as httpd:
                print(f"Serving cell analysis dashboard at http://localhost:{PORT}")
                httpd.serve_forever()

        # Start server in a separate thread
        server_thread = threading.Thread(target=start_server, daemon=True)
        server_thread.start()

        # Open the dashboard in the default browser.
        # Add the timestamp as version number to prevent loading from the cache
        webbrowser.open(f'http://localhost:{PORT}/dashboard/cell_index.html?v={time.time()}')

        # try:
        #     input("Press Enter to stop the server and close the dashboard...")
        # except KeyboardInterrupt:
        #     print("\nShutting down the server...")

    def spot_misread_density(self) -> None:
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

