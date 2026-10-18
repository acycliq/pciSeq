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
- Implements equations from the pciSeq paper
- Handles missing single-cell reference data
- Supports parallel processing via numpy operations

Dependencies:
-----------
- numpy: For numerical computations
- pandas: For data management
- scipy: For statistical operations
- numpy_groupies: For group operations
"""
import logging
from typing import Dict, List, Optional, Tuple, Union, Any

# Third-party imports
import sys
import time
import numpy as np
import numpy_groupies as npg
import pandas as pd
from scipy.special import softmax
from .utils.numba_kernels import spots_to_cell_numba_kernel
import opt_einsum as oe

# Local imports
from .datatypes.cells import Cells
from .datatypes.genes import Genes
from .datatypes.spots import Spots
from .datatypes.singleCell import SingleCell
from .datatypes.cellClass import CellClass
from .summary import collect_data
from .utils.elbo import calc_elbo
from .utils import convergence, likelihood, gaussian_model, inspection, visualisation, cell_utils
from .io import read_tsv, run_metadata
import joblib

# Configure logging
logger = logging.getLogger(__name__)


class VarBayes:
    """
    The variational Bayes model: assigns spots to cells and cells to classes.

    Parameters
    ----------
    cells_df : pd.DataFrame
        One row per cell, with columns 'label', 'area' and the centroid 'x0', 'y0',
        'z0'. stage_data produces it.
    spots_df : pd.DataFrame
        One row per spot, with columns 'x', 'y', 'z', 'plane_id', 'gene_name', 'score',
        'intensity' and 'label', the label of the cell the spot lies in, 0 for none.
        stage_data produces it.
    scRNAseq : pd.DataFrame
        Cell type definitions: mean expression per gene and cell type, genes as
        rows and cell types as columns.
    config : dict
        The resolved configuration.

    Attributes
    ----------
    cells, spots, genes, single_cell, cellTypes : object
        The parts of the model. Working with results describes the arrays on each.
    nC, nS, nG, nK : int
        Number of cells (including the background row 0), spots, genes and classes
        (including Zero).
    config : dict
        The configuration of the run: the defaults, the `opts` overrides, and the
        runtime keys `is3D`, `img_dim` and `label_map`.
    metadata : dict
        Provenance recorded when the model is built: `version`, `branch`, `commit`,
        `commit_date`, `created_at`, `python_version`, `os` and `package_versions`.
    has_converged : bool
        True when the loop stopped because the change fell below
        `CellCallTolerance`, False when it ran to `max_iter`.
    iter_delta : list of float
        The largest change in the spot assignment probabilities, one entry per
        iteration.
    iter_num : int
        Index of the last iteration run, counting from 0.
    """

    def __init__(self,
                 cells_df: pd.DataFrame,
                 spots_df: pd.DataFrame,
                 scRNAseq: pd.DataFrame,
                 config: Dict[str, Any]) -> None:
        """Initialize components and setup."""
        # Explicitly declare important instance attributes
        self.config = None
        self.iter_num = None
        self.iter_delta = []
        self.has_converged = False
        self.on_iteration_callback = None  # Optional callback for real-time visualization

        # Initialize components
        self._validate_config(config)
        self.config = config
        self._setup_components(cells_df, spots_df, scRNAseq)
        self._setup_dimensions()

        # Placeholder for other attributes
        self._scaled_exp = None
        # self._cell_explorer: Optional[CellExplorer] = None


        # Stamp this run so a result can be traced back to the code and the
        # environment that made it. Built in io/provenance.py, the one place for it.
        self.metadata = run_metadata()

    @staticmethod
    def _validate_config(config: Dict[str, Any]) -> None:
        """Check for required config parameters."""
        required = ['exclude_genes', 'max_iter', 'CellCallTolerance',
                    'rGene', 'Inefficiency', 'InsideCellBonus', 'MisreadDensity',
                    'cell_centroid_prior', 'cell_cov_prior', 'SpotReg', 'nNeighbors', 'rSpot',
                    'save_data', 'output_path', 'cell_radius', 'cell_type_prior', 'is3D',
                    ]
        missing = [param for param in required if param not in config]
        if missing:
            raise ValueError(f"Missing required config parameters: {missing}")

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
        self.cells.nbrs = self.cells.nearest_neighbours()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.genes.init_eta(self.config['rGene'], self.config['rGene'])
        # each spot's nN nearest cells, picked once from the segmentation
        # centroids. gaussian_upd moves the centroids later but this list is
        # not rebuilt, the spot keeps competing among the same cells
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)[0]
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id, self.config)
        self.cells._ini_gene_counts = np.bincount(self.spots.data.label.values, minlength=self.nC)
        self.genes._misread_density = self.genes.calc_misread_density()
        self.genes.init_rho(self.config['MisreadDensity']['default'], self._roi_volume())
        self.spots.init_gamma(self.config['rSpot'], self.config['rSpot'], [self.nC, self.nG, self.nK])
        self.init_theta()

    def _roi_volume(self) -> float:
        """Volume of the region of interest, in the same units the spots live in.

        The background misread density is N_0,g / A_total, and it gets compared
        against the cell gaussians inside the same softmax. Those gaussians are
        evaluated on anisotropy scaled coordinates, z stretched by
        voxel_size[2]/voxel_size[0], so they are densities per scaled volume.
        Counting the roi in raw voxels leaves the two sides in different units
        and the background comes out too strong by that same factor, 2.5 times
        for a 0.28/0.28/0.7 voxel, so too many spots end up as background
        instead of going to a cell.

        Only z is corrected. We take the pixels to be square in xy, which they
        are on every dataset we run. anisotropy_calc does also scale y by
        voxel_size[1]/voxel_size[0], so if a dataset ever turns up with non
        square pixels this volume will be out by that factor and it needs to go
        back in. Default voxel_size is [1, 1, 1], so 2d and isotropic runs are
        untouched.
        """
        dim = self.config['img_dim']
        vs = self.config['voxel_size']
        Sz = vs[2] / vs[0]
        return dim['w'] * dim['h'] * dim['n_planes'] * Sz

    def init_theta(self) -> None:
        geneCounts = self.cells.ini_gene_counts
        alpha = geneCounts + self.config['rTheta'] - 1



        mu = self.single_cell.mean_expression_adj + self.config['SpotReg']
        area_factor = self.cells.ini_cell_props['area_factor']
        gamma_bar = self.spots.gamma_bar
        eta_bar = self.genes.eta_bar

        beta = np.einsum('c, cgk, g, gk -> ck',
                         area_factor,
                         gamma_bar,
                         eta_bar,
                         mu) + self.config['rTheta']

        self.cells.calc_theta(alpha, beta)


    def __getstate__(self):
        """
        Get state for pickling.
        Drops the real-time viewer callback so the model can be pickled.
        """
        attributes = self.__dict__.copy()
        if 'on_iteration_callback' in attributes:
            del attributes['on_iteration_callback']
        return attributes

    @property
    def scaled_exp(self):
        """
        Get scaled expression values.

        Returns:
            np.ndarray: the scaled expression array
        """
        return self._scaled_exp


    # -------------------------------------------------------------------- #
    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        self.initialise_state()
        cell_df, gene_df = self.main_loop()
        return cell_df, gene_df

    def _step(self, name, fn):
        """Run one update step. When config['elbo_per_step'] names this step, score
        the elbo either side of it and log the change. When config['verbose'] is on,
        record its wall-time into self._step_times. When neither, this is just fn()
        with no overhead.

        Scoring here rather than on the methods themselves is deliberate. The eight
        steps carve the iteration into slices that do not overlap, so their deltas
        add up to the change over the whole iteration. Some steps call others (
        spots_to_cell runs geneCount_upd at the end), so scoring the methods would
        count those inner calls twice and the sum would not add up any more.

        The first iteration is not scored: spots.mvn_loglik_arr is None until the
        first spots_to_cell has run, and that is used here as the sign that the
        state is a fitted one and not the initialisation.
        """
        watch = self.config.get('elbo_per_step') or ()
        if (watch == 'all' or name in watch) and self.spots.mvn_loglik_arr is not None:
            before = calc_elbo(self)
            fn()
            after = calc_elbo(self)
            logger.info('elbo step %-16s %+14.2f   (%.2f -> %.2f)',
                        name, after - before, before, after)
            return
        if not self.config.get('verbose', False):
            fn()
            return
        t0 = time.perf_counter()
        fn()
        self._step_times[name] = time.perf_counter() - t0

    # -------------------------------------------------------------------- #
    def main_loop(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Run the variational Bayes loop.

        One pass updates, in this order: the gene counts per cell, the per gene
        misread density, the gene inefficiency, the cell inefficiency, gamma, the
        cell types, the spot to cell assignments, then the cell geometry. One more step, the dirichlet
        prior update, runs right after the cell types, but only when cell_type_prior
        is 'weighted'. It gives classes with more cells a higher prior; the Zero weight
        stays fixed at its cell_type_weights value. With the default, 'uniform', the
        prior stays at cell_type_weights throughout.

        The loop stops when the biggest change in the spot to cell probabilities
        drops below CellCallTolerance, or when it runs out of iterations.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: the cell dataframe with the final
            assignments and probabilities, and the gene dataframe.
        """
        p0 = None
        verbose = self.config.get('verbose', False)
        max_iter = self.config['max_iter']

        for i in range(max_iter):
            self.iter_num = i
            self._step_times = {}
            # only the verbose diagnostics read this and it is a nC by nK copy every
            # pass, so do not pay for it when nobody is going to look at it
            classProb_before = self.cells.classProb.copy() if verbose else None

            # 1. For each cell, calc the expected gene counts
            self._step('geneCount_upd', self.geneCount_upd)

            # 2. update gene-specific misread density
            self._step('rho_upd', self.rho_upd)

            # 3. calc the gene inefficiency
            self._step('eta_upd', self.eta_upd)

            # 4. calc the cell inefficiency
            self._step('theta_upd', self.theta_upd)

            # 5. calc expected gamma
            self._step('gamma_upd', self.gamma_upd)

            # 6. assign cells to cell types
            self._step('cell_to_cellType', self.cell_to_cellType)

            # 7. update the dirichlet distribution
            if self.config['cell_type_prior'] == 'weighted':
                self._step('dalpha_upd', self.dalpha_upd)

            # 8. assign spots to cells
            # spots_to_cell_numba is the fast path. spots_to_cell (the plain numpy
            # loop) is kept as the readable reference and can be swapped in on this
            # line when needed. tests/test_spots_to_cell_ab.py checks they agree.
            self._step('spots_to_cell', self.spots_to_cell_numba)

            # 9. the cell geometry: the centre and covariance of each cell's gaussian.
            # It has to come after the spot assignment, not before it. On the first
            # pass the spot weights are still the uniform 1/nN from initialise_state,
            # and the soft mean of every spot that merely has cell c among its nN
            # nearest is nowhere near cell c. Run on those, the centres jump by tens
            # of pixels and the loop never recovers. So the first spot assignment
            # uses the prior geometry, and every geometry update sees real weights.
            self._step('gaussian_upd', self.gaussian_upd)

            if verbose:
                _total = sum(self._step_times.values())
                _bd = ' '.join('%s=%.2f' % (k, v) for k, v in self._step_times.items())
                logger.info('STEP TIMES iter %d (total %.2fs): %s', i, _total, _bd)

            self.has_converged, delta = convergence.has_converged(
                self.spots, p0, self.config['CellCallTolerance'], verbose
            )
            convergence.log_iteration_diagnostics(self, i, delta, p0, classProb_before)

            # Call real-time viewer callback if provided
            if self.on_iteration_callback is not None:
                try:
                    self.on_iteration_callback(self.cells.classProb, i, delta)
                except Exception as e:
                    logger.warning(f"Real-time viewer callback failed: {e}")

            # keep track of the deltas
            self.iter_delta.append(delta)

            # replace p0 with the latest probabilities. Safe to hold the reference
            # rather than a copy because spots_to_cell always assigns a brand new
            # array, it never writes into the old one.
            p0 = self.spots.parent_cell_prob

            if self.has_converged:
                break

            if i == max_iter - 1:
                logger.info('Loop exhausted. Exiting with convergence status: %s' % self.has_converged)

        return collect_data(self.cells, self.spots, self.genes, self.config['is3D'])

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
            - self._scaled_exp: scaled expression
            - self.spots._log_gamma_bar: Log of expected gamma values
            - self.spots._gamma_bar: Expected gamma values
        """
        cells = self.cells
        cfg = self.config

        self._scaled_exp = likelihood.scaled_exp(cells.ini_cell_props['area_factor'],
                                            self.single_cell.mean_expression_adj.values)

        beta = self.scaled_exp * self.genes.eta_bar[:, None] * self.cells.theta_bar[:,None, :]+ cfg['rSpot']
        rho = cfg['rSpot'] + cells.geneCount

        self.spots._post_shape = rho
        self.spots._post_rate = beta
        self.spots._log_gamma_bar = self.spots.logGammaExpectation(rho, beta)
        self.spots._gamma_bar = self.spots.gammaExpectation(rho, beta)

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
        contr = likelihood.compute_gene_loglikelihood_matrix(self)

        # per cell, gene and class contributions to the NB loglik from this class update.
        # check_cell reads it, so dont remove it.
        self.cells.nb_contr = contr
        contr = np.sum(contr, axis=1)
        mrf = self.cells.calc_mrf()
        # mrf = self.cells.classProb[self.cells.nbrs].sum(axis=1)
        # keep the mrf around for the diagnostics export
        self.cells.mrf = mrf
        wCellClass = contr + self.cellTypes.log_prior + mrf
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

        # Gene-specific misread density (learned per gene)
        log_rho = self.genes.log_rho_bar[self.spots.gene_id]

        # pre-populate last column
        wSpotCell[:, -1] = log_rho
        mvn_loglik_arr = np.zeros(wSpotCell.shape)
        attention = np.zeros(wSpotCell.shape)
        expr_fluctuations = np.zeros(wSpotCell.shape)
        cell_inefficiency = np.zeros(wSpotCell.shape)
        gene_inefficiency = np.zeros(wSpotCell.shape)

        # Materialize once before the loop (same for all neighbors)
        log_gamma_bar_arr = self.spots.log_gamma_bar

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]
            log_theta_bar = np.log(self.cells.theta_bar[sn])

            # multiply and sum over cells. In practice this means that when high expected counts
            # are aligned with high cell class probs this term will be high
            term_1 = np.einsum('ij, ij -> i', expected_counts, cp)

            log_gamma_bar = log_gamma_bar_arr[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)

            term_3 = np.einsum('ij, ij -> i', cp, log_theta_bar)

            # wSpotCell[:, n] = term_1 + term_2 + logeta_bar + loglik[:, n]
            mvn_loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells, self.config['is3D'])
            wSpotCell[:, n] = term_1 + term_2 + term_3 + logeta_bar + mvn_loglik
            mvn_loglik_arr[:, n] = mvn_loglik
            attention[:, n] = term_1
            expr_fluctuations[:, n] = term_2
            cell_inefficiency[:, n] = term_3
            gene_inefficiency[:, n] = logeta_bar

        # apply inside cell bonus
        bonus_mask = self.spots.bonus_mask * self.config['InsideCellBonus']
        wSpotCell += bonus_mask

        # update the prob a spot belongs to a neighboring cell
        self.spots.parent_cell_prob = softmax(wSpotCell, axis=1)
        self.spots.mvn_loglik_arr = mvn_loglik_arr
        self.spots.attention = attention
        self.spots.expr_fluctuations = expr_fluctuations
        self.spots.cell_inefficiency = cell_inefficiency
        self.spots.gene_inefficiency = gene_inefficiency

        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()

    # -------------------------------------------------------------------- #
    def spots_to_cell_numba(self) -> None:
        """Compute the spot-to-cell assignment probabilities (numba version).

        This is the fast path used by default. spots_to_cell is the equivalent plain
        numpy loop, kept as the readable reference; the two must give the same result,
        which is checked in tests/test_spots_to_cell_ab.py.

        Returns
        -------
        None
            Sets self.spots.parent_cell_prob and the per-spot diagnostic arrays
            (attention, expr_fluctuations, cell_inefficiency, gene_inefficiency,
            mvn_loglik_arr).

        Notes
        -----
        The three score terms are computed by spots_to_cell_numba_kernel in
        numba_kernels.mvn_loglik is computed the same way as in spots_to_cell.
        """
        nN = self.nN
        nNb = nN - 1
        nS = self.spots.data.gene_name.shape[0]

        gn = self.spots.data.gene_name.values
        expected_counts = np.ascontiguousarray(self.single_cell.log_mean_expression.loc[gn].values)  # [nS, nK]
        logeta_bar = np.ascontiguousarray(self.genes.logeta_bar[self.spots.gene_id])     # [nS]
        log_rho = self.genes.log_rho_bar[self.spots.gene_id]                             # [nS]
        log_gamma_bar_arr = np.ascontiguousarray(self.spots.log_gamma_bar)     # [nC, nG, nK]
        log_theta_bar_all = np.ascontiguousarray(np.log(self.cells.theta_bar))           # [nC, nK]
        classProb = np.ascontiguousarray(self.cells.classProb)                           # [nC, nK]
        parent = np.ascontiguousarray(self.spots.parent_cell_id)                         # [nS, nN]
        gene_id = np.ascontiguousarray(self.spots.gene_id)                               # [nS]

        wSpotCell = np.zeros([nS, nN], dtype=np.float64)
        wSpotCell[:, -1] = log_rho
        mvn_loglik_arr = np.zeros([nS, nN])
        attention = np.zeros([nS, nN])
        expr_fluctuations = np.zeros([nS, nN])
        cell_inefficiency = np.zeros([nS, nN])
        gene_inefficiency = np.zeros([nS, nN])

        # numba kernel fills in term_1/2/3 for all spots (see numba_kernels)
        spots_to_cell_numba_kernel(parent, gene_id, classProb, log_gamma_bar_arr, log_theta_bar_all,
                                   expected_counts, logeta_bar, nNb,
                                   wSpotCell, attention, expr_fluctuations, cell_inefficiency, gene_inefficiency)

        # mvn_loglik is not in the kernel yet; compute per neighbour and add it in
        for n in range(nNb):
            sn = parent[:, n]
            mvn = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells, self.config['is3D'])
            wSpotCell[:, n] += mvn
            mvn_loglik_arr[:, n] = mvn

        bonus_mask = self.spots.bonus_mask * self.config['InsideCellBonus']
        wSpotCell += bonus_mask

        self.spots.parent_cell_prob = softmax(wSpotCell, axis=1)
        self.spots.mvn_loglik_arr = mvn_loglik_arr
        self.spots.attention = attention
        self.spots.expr_fluctuations = expr_fluctuations
        self.spots.cell_inefficiency = cell_inefficiency
        self.spots.gene_inefficiency = gene_inefficiency

        # Since the spot-to-cell assignments changed you need to update the gene counts now.
        # Same as the plain spots_to_cell above. Keep it, it is not redundant.
        self.geneCount_upd()


    # -------------------------------------------------------------------- #
    def rho_upd(self) -> None:
        """Updates gene-specific misread density (rho_g).

        Uses the expected number of background spots per gene to update
        the Gamma posterior for each gene's misread density.
        """
        background_counts = np.bincount(
            self.spots.gene_id,
            self.spots.parent_cell_prob[:, -1],
            minlength=self.nG
        )
        self.genes.calc_rho(background_counts)
        if self.config.get('verbose', False):
            logger.info(f"rho_upd: bg_counts min/max={background_counts.min():.1f}/{background_counts.max():.1f}, "
                         f"rho_bar min/max={self.genes.rho_bar.min():.2e}/{self.genes.rho_bar.max():.2e}, "
                         f"log_rho min/max={self.genes.log_rho_bar.min():.4f}/{self.genes.log_rho_bar.max():.4f}")

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
        gamma_bar = self.spots.gamma_bar
        theta_bar = self.cells.theta_bar

        # same thing as self.spots.zero_class_counts(self.spots.gene_id, classProb[:, -1]),
        # checked on espio 14-Sep-2026, 40 iterations, same deltas to float32 noise
        zero_class_counts = oe.contract('c, cg -> g', classProb[:, -1], self.cells.geneCount, optimize='optimal')

        # Calcs the sum in the Gamma distribution (equation 5). The zero class
        # is excluded from the sum, hence the arrays in the einsum below stop at :-1
        # Note. Row 0 is the placeholder for the background, not a cell, and strictly it
        # should not be in the sum. It is left in. geneCount_upd zeroes its counts, so the
        # model calls it Zero almost entirely and what it adds here is next to nothing.
        class_total_counts = oe.contract('ck, gk, c, cgk, ck -> g',
                                         classProb[:, :-1],
                                         mu.values[:, :-1],
                                         area_factor,
                                         gamma_bar[:, :, :-1],
                                         theta_bar[:,:-1], optimize='optimal')
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
        """Update the cell geometry, the NIW posterior over each cell's centroid
        and covariance. The maths is on the-model/gaussian.md.

        centroid_upd does kappa_hat and m_hat, cov_upd does nu_hat and V_hat.
        Both need the soft count, mean and scatter of the cell's spots so those
        are worked out once here and passed down.

        Two things are left as they were, on purpose. The nN nearest cells of
        each spot were picked from the segmentation centroids in
        initialise_state and are not picked again when the centres move. And
        the inside cell bonus keeps its area factor from read_image_objects,
        the one for a round cell of the mean radius; that the gaussian may now
        be stretched or tilted is ignored there.
        """
        n_bar, x_bar, S = gaussian_model.soft_stats(self.spots, self.cells, self.cells.dim)
        self.centroid_upd(n_bar, x_bar)
        self.cov_upd(n_bar, x_bar, S)

    # -------------------------------------------------------------------- #
    def centroid_upd(self, n_bar, x_bar) -> None:
        """kappa_hat = kappa0 + n_bar and m_hat, the weighted average of the
        segmentation centroid (worth kappa0 spots) and the soft mean of the spots
        (worth n_bar). Cell 0 is the background, it has n_bar = 0 and stays put.
        In 2D only x and y move, z stays whatever the segmentation said.
        """
        kappa0 = self.config['cell_centroid_prior']['default']
        d = self.cells.dim
        m0 = self.cells.ini_centroids().values[:, :d].astype(np.float64)

        kappa_hat = kappa0 + n_bar
        m_hat = (kappa0 * m0 + n_bar[:, None] * x_bar) / kappa_hat[:, None]

        centroid = self.cells.ini_centroids()
        centroid.iloc[:, :d] = m_hat.astype(np.float32)
        self.cells.centroid = centroid
        self.cells.kappa = kappa_hat.astype(np.float32)

    # -------------------------------------------------------------------- #
    def cov_upd(self, n_bar, x_bar, S) -> None:
        """nu_hat = nu0 + n_bar and the scale

            V_hat = V0 + S + kappa0 n_bar / (kappa0 + n_bar) (x_bar - m0)(x_bar - m0)^T

        S is the scatter about the soft mean x_bar, the last term is what the
        centre having moved away from m0 adds (the parallel axis theorem, see
        the docs page). cov gets E[Sigma] = V_hat / (nu_hat - d - 1) for the
        ellipses, and the eigen pairs of V_hat are what the position term uses.
        """
        kappa0 = self.config['cell_centroid_prior']['default']
        nu0 = self.config['cell_cov_prior']['default']
        d = self.cells.dim
        m0 = self.cells.ini_centroids().values[:, :d].astype(np.float64)

        dm = x_bar - m0
        shift = (kappa0 * n_bar / (kappa0 + n_bar))[:, None, None] * oe.contract('ci,cj->cij', dm, dm)
        V_hat = self.cells.ini_scale() + S + shift
        nu_hat = nu0 + n_bar

        self.cells.scale = V_hat.astype(np.float32)
        self.cells.nu = nu_hat.astype(np.float32)
        self.cells.cov = (V_hat / (nu_hat - d - 1)[:, None, None]).astype(np.float32)
        self.cells.eig_vals, self.cells.eig_vecs = np.linalg.eigh(self.cells.scale)

    # -------------------------------------------------------------------- #
    def dalpha_upd(self) -> None:
        """
        Updates cell type prior distribution parameters.

        Adjusts Dirichlet parameters based on:
            1. Current cell type assignments
            2. Initial alpha values
        """
        # logger.info('Update cell type (marginal) distribution')
        # Only update real classes (exclude Zero, which is the last column)
        zeta_real = self.cells.classProb[:, :-1].sum(axis=0)
        alpha_0 = np.ones(self.cellTypes.nK - 1, dtype=np.float32)
        self.cellTypes.alpha = zeta_real + alpha_0

    # -------------------------------------------------------------------- #
    def theta_upd(self):
        geneCounts = self.cells.geneCount.sum(axis=1)
        alpha = geneCounts + self.config['rTheta'] - 1

        mu = self.single_cell.mean_expression_adj + self.config['SpotReg']
        area_factor = self.cells.ini_cell_props['area_factor']
        gamma_bar = self.spots.gamma_bar
        eta_bar = self.genes.eta_bar

        beta = np.einsum('c, cgk, g, gk -> ck',
                         area_factor,
                         gamma_bar,
                         eta_bar,
                         mu) + self.config['rTheta']

        self.cells.calc_theta(alpha, beta)

    # -------------------------------------------------------------------- #
    def heatmap_counts_per_class(self):
        """Display the interactive heatmap."""
        return visualisation.heatmap_counts_per_class(self)

    def calculate_genes_log_likelihood_contr(self, label):
        return likelihood.calculate_genes_log_likelihood_contr(self, label)

    def check_cell(self, my_label, user_class, top_n=10, show_plot=True, top_classes=5):
        """
        Compare the assigned cell type of a cell with another cell type.

        The assigned type is the type with the highest probability in classProb. The
        per-gene log-likelihoods, log prior and spatial (MRF) term are those used in
        the last cell type update. The figure shows the genes that most favour each
        type, the three score components for both types, and the posterior over all
        types for the most likely ones.

        Parameters
        ----------
        my_label : int
            Cell label, as in the input segmentation.
        user_class : str
            Cell type to compare against. Must differ from the assigned type.
        top_n : int, default 10
            Maximum number of genes shown on each side: the genes that most favour the
            assigned type and those that most favour user_class. Only genes with a
            nonzero difference in that direction are shown, so a side can have fewer.
        show_plot : bool, default True
            Draw the figure.
        top_classes : int, default 5
            Number of types shown in the posterior chart, ordered by probability.
            user_class is added if it is not among them.

        Returns
        -------
        gene_expression_data : pd.DataFrame
            One row per selected gene. Columns are the mean counts of the gene in cells
            of each type, the expected count of the gene in this cell under each type,
            and its counts in this cell. The counts are sums of spot assignment
            probabilities, so they are estimates, not measurements.
        contr : pd.DataFrame
            Log-likelihood of each gene in this cell under the two types, and their
            difference, 'diff'.
        fig : matplotlib.figure.Figure or None
            The figure, or None if show_plot is False.

        Raises
        ------
        ValueError
            If user_class is not a cell type, or is the assigned type of the cell.
        """
        return inspection.check_cell(self, my_label, user_class, top_n, show_plot, top_classes)

    def check_spot(self, spot_id, show_plot=True):
        """
        Break down the assignment of a spot to its candidate cells.

        The score of each candidate cell is the sum of the terms from the last
        spot-to-cell update: the spatial log-likelihood ('Gaussian fit'), the expected log mean
        expression under the cell's class probabilities ('class expression'), the expected
        log theta ('cell scale'), the expected log gamma ('cell-gene scale'), the log eta
        of the spot's gene ('gene efficiency') and the inside-cell bonus ('bonus'). The
        score of the background is the log misread density of the gene. Probabilities are
        the softmax of the scores.

        Parameters
        ----------
        spot_id : int
            Spot id, the index of the spots table.
        show_plot : bool, default True
            Draw the score and probability charts.

        Returns
        -------
        pd.DataFrame
            One row per candidate cell, indexed by cell label, and a final 'background'
            row. Columns are the score terms above,
            'misread', 'sum' and 'prob'.
        """
        return inspection.check_spot(self, spot_id, show_plot)

    def to_internal(self, label):
        """
        Convert a segmentation cell label to the internal one.

        The internal label is the cell's row index in the arrays in the pickle file,
        1 to nC - 1, with row 0 the background. It differs from the segmentation label only when the
        input labels were not sequential, which usually happens because `remove_flat_cells`
        dropped the cells that span a single plane.

        Parameters
        ----------
        label : int or list of int
            Cell label, as in the input segmentation.

        Returns
        -------
        int or list of int
            The row to index the arrays with. The same value back when no renumbering
            took place.

        Raises
        ------
        KeyError
            If the label is not in the segmentation.
        """
        return cell_utils.to_internal(label, self.config['label_map'])

    def to_external(self, label):
        """
        Convert an internal cell label to the segmentation one.

        The opposite of `to_internal`. Use it on a row index taken from an array in the
        pickle file or from diagnostics.db, to get the label the segmentation gave that cell.

        Parameters
        ----------
        label : int or list of int
            Row index into the arrays in the pickle file.

        Returns
        -------
        int or list of int
            The segmentation label. The same value back when no renumbering took place.

        Raises
        ------
        KeyError
            If the row index is out of range.
        """
        return cell_utils.to_external(label, self.config['label_map'])

    def read_tsv(self, filepath):
        """
        Read a tsv file written by pciSeq into a DataFrame.

        Columns that hold lists or dicts are parsed back from text.

        Parameters
        ----------
        filepath : str
            Path to the file, e.g. cellData.tsv or geneData.tsv.

        Returns
        -------
        pd.DataFrame
        """
        return read_tsv(filepath)

    # def trellis_plot(self, label, flatfile_folder):
    #     return visualisation.trellis_plot(self, label, flatfile_folder)
