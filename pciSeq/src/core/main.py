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

import time
import logging
from typing import Dict, List, Optional, Tuple, Union

# Third-party imports
import numpy as np
import numpy_groupies as npg
import pandas as pd
from dask.delayed import delayed
import scipy.spatial as spatial
from scipy.special import softmax
from typing import Dict, Any

# Local imports
from .datatypes import Cells, Spots, Genes, SingleCell, CellType
from .summary import collect_data
from .analysis import CellExplorer
from . import utils
from ...src.diagnostics.redis_publisher import RedisPublisher
from ...src.diagnostics.utils import RedisDB

main_logger = logging.getLogger(__name__)


class VarBayes:
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

        # Placeholder for other attributes
        self._scaled_exp = None
        self._cell_explorer: Optional[CellExplorer] = None

    def _validate_config(self, config: Dict[str, Any]) -> None:
        """Check for required config parameters."""
        required = ['exclude_genes', 'max_iter', 'CellCallTolerance',
                    'rGene', 'Inefficiency', 'InsideCellBonus', 'MisreadDensity',
                    'SpotReg', 'nNeighbors', 'rSpot', 'save_data', 'output_path',
                    'launch_viewer', 'launch_diagnostics', 'is_redis_running',
                    'cell_radius', 'cell_type_prior', 'mean_gene_counts_per_class',
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
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)
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
        p0 = None
        cell_df = None
        gene_df = None
        max_iter = self.config['max_iter']

        self.initialise_state()
        for i in range(max_iter):
            self.iter_num = i

            # 1. For each cell, calc the expected gene counts
            self.geneCount_upd()

            # 2. calc expected gamma
            self.gamma_upd()

            # 3 update correlation matrix and variance of the gaussian distribution
            if self.single_cell.isMissing or (self.config['InsideCellBonus'] is False):
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

            self.has_converged, delta = utils.has_converged(self.spots, p0, self.config['CellCallTolerance'])
            main_logger.info('Iteration %d, mean prob change %f' % (i, delta))

            # keep track of the deltas
            self.iter_delta.append(delta)

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            if self.redis_publisher:
                self.redis_publisher.publish_diagnostics(self, self.iter_num, self.has_converged)

            if self.has_converged:
                # self.counts_within_radius(20)
                self.cell_explorer.view_cell(2259)
                cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell)
                break

            if i == max_iter - 1:
                main_logger.info('Loop exhausted. Exiting with convergence status: %s' % self.has_converged)
        return cell_df, gene_df

    # -------------------------------------------------------------------- #
    def geneCount_upd(self):
        """
        Produces a matrix numCells-by-numGenes where element at position (c,g) keeps the expected
        counts of gene g  in cell c.
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
    def gamma_upd(self):
        """
         Implements equation (3) of the Qian paper
        """
        cells = self.cells
        cfg = self.config

        self._scaled_exp = delayed(utils.scaled_exp(cells.ini_cell_props['area_factor'],
                                            self.single_cell.mean_expression.values,
                                            self.genes.eta_bar))

        beta = self.scaled_exp.compute() + cfg['rSpot']
        rho = cfg['rSpot'] + cells.geneCount

        self.spots._gamma_bar = self.spots.gammaExpectation(rho, beta)
        self.spots._log_gamma_bar = self.spots.logGammaExpectation(rho, beta)

    # -------------------------------------------------------------------- #
    def cell_to_cellType(self):
        """
        returns an array of size numCells-by-numCellTypes where element in position [i,j]
        keeps the probability that cell i has cell type j
        Implements equation (2) of the Qian paper
        :param spots:
        :param config:
        :return:
        """

        ScaledExp = self.scaled_exp.compute()
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negative_binomial_loglikelihood(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cellTypes.log_prior
        pCellClass = softmax(wCellClass, axis=1)
        del contr

        self.cells.classProb = pCellClass

    # -------------------------------------------------------------------- #
    def spots_to_cell(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        nN = self.nN
        nS = self.spots.data.gene_name.shape[0]

        wSpotCell = np.zeros([nS, nN], dtype=np.float64)
        gn = self.spots.data.gene_name.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values
        logeta_bar = self.genes.logeta_bar[self.spots.gene_id]

        # pre-populate last column
        misread = self.spot_misread_density()
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
            mvn_loglik = self.spots.mvn_loglik(self.spots.xy_coords, sn, self.cells)
            wSpotCell[:, n] = term_1 + term_2 + logeta_bar + mvn_loglik

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
    def eta_upd(self):
        """
        Calcs the expected eta
        Implements equation (5) of the Qian paper
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
    def gaussian_upd(self):
        self.centroid_upd()
        self.cov_upd()

    # -------------------------------------------------------------------- #
    def centroid_upd(self):
        spots = self.spots

        # get the total gene counts per cell
        N_c = self.cells.total_counts

        xy_spots = spots.xy_coords
        prob = spots.parent_cell_prob
        n = self.cells.config['nNeighbors'] + 1

        # multiply the x coord of the spots by the cell prob
        a = np.tile(xy_spots[:, 0], (n, 1)).T * prob

        # multiply the y coord of the spots by the cell prob
        b = np.tile(xy_spots[:, 1], (n, 1)).T * prob

        # aggregated x and y coordinate
        idx = spots.parent_cell_id
        x_agg = npg.aggregate(idx.ravel(), a.ravel(), size=len(N_c))
        y_agg = npg.aggregate(idx.ravel(), b.ravel(), size=len(N_c))

        # get the estimated cell centers
        x_bar = np.nan * np.ones(N_c.shape)
        y_bar = np.nan * np.ones(N_c.shape)

        x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
        y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]

        # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
        xy_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T)))

        # if you have a value for the estimated centroid use that, otherwise
        # use the initial (starting values) centroids
        ini_cent = self.cells.ini_centroids()
        xy_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y']])))

        # # sanity check. NaNs or Infs should appear together
        # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
        # use the fitted centroids where possible otherwise use the initial ones
        xy_bar[np.isfinite(x_bar)] = xy_bar_fitted[np.isfinite(x_bar)]
        self.cells.centroid = pd.DataFrame(xy_bar, columns=['x', 'y'])
        # print(np.array(list(zip(x_bar.T, y_bar.T))))

    # -------------------------------------------------------------------- #
    def cov_upd(self):
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
    def mu_upd(self):
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
    def dalpha_upd(self):
        # main_logger.info('Update cell type (marginal) distribution')
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
    def spot_misread_density(self):
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


