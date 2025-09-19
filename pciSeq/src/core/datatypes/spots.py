# Standard library imports
import logging
from typing import Tuple, Dict, Any

# Third party imports
import numpy as np
import pandas as pd
import scipy
import opt_einsum as oe

from ..utils.ops_utils import gene_density

spots_logger = logging.getLogger(__name__)


class Spots(object):
    """
    Processes RNA spot detection data, including spatial coordinates and
    cell assignments. Includes methods for reading spot data, calculating
    nearest cells, and managing gene counts.

    Attributes:
        config (dict): Configuration parameters for spot data.
        data (pd.DataFrame): DataFrame containing spot data.
        data_excluded (pd.DataFrame): DataFrame containing excluded spot data.
        nS (int): Number of spots.
        unique_gene_names (np.array): Array of unique gene names.
        _gamma_bar (np.array): Gamma bar values for spots.
        _log_gamma_bar (np.array): Log gamma bar values for spots.
        _gene_id (np.array): Gene IDs for spots.
        _counts_per_gene (np.array): Counts per gene for spots.
    """

    def __init__(self, spots_df: pd.DataFrame, config: Dict):
        """
        Initializes the Spots object with spot data and configuration.

        Parameters:
            data (pd.DataFrame): DataFrame containing spot data.
            config (dict): Configuration parameters for spot data.
        """
        self._parent_cell_prob = None
        self._parent_cell_id = None
        self.Dist = None
        self.config = config
        self.data, self.data_excluded = self.read(spots_df)
        self.nS = self.data.shape[0]
        self.unique_gene_names = None
        self._gamma_bar = None
        self._log_gamma_bar = None
        self._gene_id = None
        self._counts_per_gene = None
        [_, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True,
                                                            return_counts=True)
        self._plane_adj = None
        self.mvn_loglik_arr = None
        self.attention = None
        self.expr_fluctuations = None

    def __getstate__(self):
        """
        Customizes the state for pickling, excluding certain attributes.

        Returns:
            dict: Attributes to be serialized.
        """
        attributes = self.__dict__.copy()
        del attributes['_gamma_bar']
        del attributes['_log_gamma_bar']
        return attributes

    # ---------------- PROPERTIES ---------------- #
    @property
    def gene_id(self) -> np.ndarray:
        """Returns the gene IDs for spots."""
        return self._gene_id

    @gene_id.setter
    def gene_id(self, val: np.ndarray):
        """Sets the gene IDs for spots."""
        self._gene_id = val.astype(np.int32)

    @property
    def counts_per_gene(self) -> np.ndarray:
        """Returns the counts per gene for spots."""
        return self._counts_per_gene

    @counts_per_gene.setter
    def counts_per_gene(self, val: np.ndarray):
        """Sets the counts per gene for spots."""
        self._counts_per_gene = val.astype(np.int32)

    @property
    def gamma_bar(self) -> np.ndarray:
        """Returns the gamma bar values for spots."""
        return self._gamma_bar

    @property
    def log_gamma_bar(self) -> np.ndarray:
        """Returns the log gamma bar values for spots."""
        return self._log_gamma_bar

    @property
    def xyz_coords(self) -> np.ndarray:
        """Returns the spatial coordinates of spots."""
        return self.data[['x', 'y', 'z']].values.astype(np.float32)

    @property
    def parent_cell_prob(self):
        return self._parent_cell_prob

    @parent_cell_prob.setter
    def parent_cell_prob(self, val: np.ndarray):
        """Sets the parent cell probabilities for spots."""
        self._parent_cell_prob = val

    @property
    def parent_cell_id(self) -> np.ndarray:
        """Returns the parent cell IDs for spots."""
        return self._parent_cell_id

    @parent_cell_id.setter
    def parent_cell_id(self, val: np.ndarray):
        """Sets the parent cell IDs for spots."""
        self._parent_cell_id = val.astype(np.uint32)

    @property
    def bonus_mask(self):
        """
        Compute a boolean mask for bonus assignment in spot-to-cell matching.

        This property returns a boolean array of shape (n_spots, n_neighbors) where:
          - Each row corresponds to a spot.
          - Each column corresponds to a neighboring cell.
          - A value of True at [i, j] indicates that the label of spot i (from self.data.label)
            matches the parent's cell ID for neighboring cell j (from self.parent_cell_id). This
            implies that spot i lies within the boundaries of cell j and should receive a bonus in
            the spot-to-cell assignment.
          - The last column represents the background and is always set to False (i.e., no bonus).

        Returns:
            np.ndarray: Boolean mask of shape (n_spots, n_neighbors).
        """
        # Get the spot labels as a column vector (shape: n_spots x 1)
        cell_label = self.data.label.values[:, None]

        # Get the parent's cell IDs (shape: n_spots x n_neighbors)
        parent_cell_id = self.parent_cell_id

        # Use broadcasting to compare the spot labels against each parent's cell ID.
        bonus_mask = (cell_label == parent_cell_id)

        # Set the last column (background) to False (no bonus).
        bonus_mask[:, -1] = False

        return bonus_mask

    @property
    def plane_adj(self):
        return self._plane_adj

    # ---------------- METHODS ---------------- #
    def read(self, spots_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Reads and processes spot data, excluding specified genes.

        Parameters:
            spots_df (pd.DataFrame): DataFrame containing spot data.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Processed spot data and excluded spot data.
        """
        # No need for x_global, y_global to be in the spots_df at first place.
        # Instead of renaming here, you could just use 'x' and 'y' when you
        # created the spots_df
        # spots_df = spots_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

        # remove a gene if it is on the exclude list
        exclude_genes = self.config['exclude_genes']
        gene_mask = [True if d not in exclude_genes else False for d in spots_df.gene_name]
        neg_gene_mask = [True if d in exclude_genes else False for d in spots_df.gene_name]
        spots_copy = spots_df.copy()
        spots_df = spots_copy.loc[gene_mask]
        spots_excluded_df = spots_copy.loc[neg_gene_mask]
        return spots_df, spots_excluded_df

    def cells_nearby(self, cells: 'Cells') -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds nearby cells for each spot.

        Parameters:
            cells (Cells): Cells object containing cell data.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Neighbors and cell probabilities.
        """
        spotZYX = self.data[['z', 'y', 'x']].values

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        Dist, neighbors = nbrs.kneighbors(spotZYX)
        self.Dist = Dist.astype(np.float32)

        # last column is for misreads.
        neighbors[:, -1] = 0

        # make an array assigning 100% prob of any given cell belonging to its closest neighbour
        cellProb = np.zeros(neighbors.shape, dtype=np.float32)
        cellProb[:, 0] = np.ones(neighbors.shape[0])

        # the second return value is not getting used. maybe in the future
        return neighbors.astype(np.int32), cellProb

    def ini_cellProb(self, neighbors: np.ndarray, cfg: Dict[str, Any]) -> np.ndarray:
        nS = self.data.shape[0]
        nN = cfg['nNeighbors'] + 1
        SpotInCell = self.data.label
        # assert (np.all(SpotInCell.index == neighbors.index))

        ## sanity check (this actually needs to be rewritten)
        # mask = np.greater(SpotInCell, 0, where=~np.isnan(SpotInCell))
        # sanity_check = neighbors[mask, 0] + 1 == SpotInCell[mask]
        # assert ~any(sanity_check), "a spot is in a cell not closest neighbor!"

        pSpotNeighb = np.zeros([nS, nN], dtype=np.float32)
        pSpotNeighb[neighbors == SpotInCell.values[:, None]] = 1
        pSpotNeighb[SpotInCell == 0, -1] = 1

        # you might have case where the sum across all columns is Zero.
        # That can happen if for example the spot is inside the cell boundaries
        # of a cell, that cell however is not one of nN-th closest cells.
        # For example a cell spans the full 3d stack, hence the centroid would be
        # around the mid-plane somewhere. You have however a spot on one of the
        # first few planes and lets say in it is inside this big cell that spans the full z-stack.
        # Assume also that we have some cells which span only some of the first few planes of
        # the stack. Their centroids could be closer to the spot than the centroid of the cell the
        # spot lies within.
        # In these case assign the spot to the background
        mask = pSpotNeighb.sum(axis=1)
        pSpotNeighb[mask == 0, -1] = 1

        ## Add a couple of checks here
        return pSpotNeighb

    # def loglik(self, cells, cfg):
    #     # area = cells.ini_cell_props['area'][1:]
    #     # mcr = np.mean(np.sqrt(area / np.pi)) * 0.5  # This is the meanCellRadius
    #     mcr = cells.mcr
    #     dim = 2  # dimensions of the normal distribution: Bivariate
    #     # Assume a bivariate normal and calc the likelihood
    #     D = -self.Dist ** 2 / (2 * mcr ** 2) - dim/2 * np.log(2 * np.pi * mcr ** 2)
    #
    #     # last column (nN-closest) keeps the misreads,
    #     D[:, -1] = np.log(cfg['MisreadDensity'])
    #
    #     mask = np.greater(self.data.label.values, 0, where=~np.isnan(self.data.label.values))
    #     D[mask, 0] = D[mask, 0] + cfg['InsideCellBonus']
    #     return D

    def mvn_loglik(self, data, cell_label, cells, is3D):
        """
        Calculates the multivariate normal log likelihood for spots.

        Parameters:
            data (np.array): Spot data.
            cell_label (np.array): Cell labels for spots.
            cells (Cells): Cells object containing cell data.
            is3D (bool): Whether the data is 3D.

        Returns:
            np.array: Log likelihood values.
        """
        centroids = cells.centroid.values[cell_label]
        covs = cells.cov[cell_label]
        eig_vals = cells.eig_vals[cell_label]
        eig_vecs = cells.eig_vecs[cell_label]
        if not is3D:
            data = data[:, :-1]
            centroids = centroids[:, :-1]
        out = self.multiple_logpdfs(data, centroids, covs, eig_vals, eig_vecs)
        return out

    def multiple_logpdfs(self, x: np.ndarray, means: np.ndarray, cov: np.ndarray, vals: np.ndarray, vecs: np.ndarray) -> np.ndarray:
        """
        vectorised mvn log likelihood evaluated at multiple pairs of (centroid_1, cov_1), ..., (centroid_N, cov_N)
        Taken from http://gregorygundersen.com/blog/2020/12/12/group-multivariate-normal-pdf/
        """
        # Thankfully, NumPy broadcasts `eigh`.
        # vals, vecs = np.linalg.eigh(covs)

        # Compute the log determinants across the second axis.
        logdets = np.sum(np.log(vals), axis=1)

        # Invert the eigenvalues.
        valsinvs = 1. / vals

        # Add a dimension to `valsinvs` so that NumPy broadcasts appropriately.
        Us = vecs * np.sqrt(valsinvs)[:, None]
        devs = x - means

        # Use `einsum` for matrix-vector multiplications across the first dimension.
        devUs = oe.contract('ni,nij->nj', devs, Us, optimize='optimal')

        # Compute the Mahalanobis distance by squaring each term and summing.
        mahas = np.sum(np.square(devUs), axis=1)

        # Compute and broadcast scalar normalizers.
        dim = len(vals[0])
        log2pi = np.log(2 * np.pi)

        return -0.5 * (dim * log2pi + mahas + logdets)

    def zero_class_counts(self, geneNo, pCellZero):
        """
        *****************************
         ******** DEPRECATED ++++++++
         ****** TO BE REMOVED *******
         ****************************
        Calculates gene counts for the zero expressing class.

        Parameters:
            geneNo (np.array): Gene numbers for spots.
            pCellZero (np.array): Probabilities of zero expression for cells.

        Returns:
            np.array: Total predicted zero counts per gene.
        """
        # for each spot get the ids of the 3 nearest cells
        spotNeighbours = self.parent_cell_id[:, :-1]

        # get the corresponding probabilities
        neighbourProb = self.parent_cell_prob[:, :-1]

        # prob that a spot belongs to a zero expressing cell
        pSpotZero = np.sum(neighbourProb * pCellZero[spotNeighbours], axis=1)

        # aggregate per gene id
        TotPredictedZ = np.bincount(geneNo, pSpotZero)
        return TotPredictedZ

    def gammaExpectation(self, rho, beta):
        """
        Calculates the expectation of a gamma distribution

        In this context:
          - rho represents gene counts for each cell and has shape (nC, nG),
            where nG is the number of genes and nC is the number of cells.
          - beta represents (scaled) gene expression from scRNAseq and has shape (nG, nK)

        The expectation is computed so that for each cell c, gene g, and class k:
            result[c, g, k] = rho[c, g] / beta[g, k]

        This is achieved using an einsum operation that is equivalent to:
            rho[:, :, None] / beta
        as verified by:
            np.allclose(rho[:, :, None] / beta, np.einsum('cg, gk -> cgk', rho, 1 / beta))

        Parameters:
            rho (np.array): Gene counts per cell with shape (nG, nC).
            beta (np.array): Scaled gene expression values with shape (nC, nG, nK) or (nG, nK)

        Returns:
            np.array: Expected gamma values with shape (nC, nG, nK).
        """

        if len(beta.shape) == 3:
            subscripts = 'cg, cgk -> cgk'
        else:
            subscripts = 'cg, gk -> cgk'
        return oe.contract(subscripts, rho, 1 / beta,  optimize='optimal')

    def logGammaExpectation(self, rho, beta):
        """
        Calculates the log expectation of a gamma distribution.

        Parameters:
            rho (np.array): Shape parameters.
            beta (np.array): Rate parameters.

        Returns:
            np.array: Log expected values.
        """
        r = rho[:, :, None]
        return scipy.special.psi(r) - np.log(beta)

    def misread_density(self, genes):
        """
        convenience functon that takes the genes misread density and
        aligns to the spots. Every spot will have a misread value based
        on its gene
        """
        misread_dict = genes.misread_density.to_dict()

        # Convert to array and align directly with spots
        v = np.array(list(misread_dict.values()))
        v = v[self.gene_id]  # Align with spots
        return v


    def calc_plane_adj(self, cells, config):
        """
        Returns an array of shape (nS, nN) where each row corresponds to the
        adjustment factor that needs to be applied to the single cell data for
        given spot and neighboring cell.
        The experimental data are not uniformly distributed across the z-stack.
        You might have cases where the lower planes contain very few gene reads
        but the higher planes contain a lot of gene reads for the same gene.
        Hence we need to adjust the single cell data for each gene depending on
        the plane and the likely parent cell it belongs to.
        Parameters
        ----------
        cells
        config

        Returns
        -------

        """
        density = gene_density(self, config)

        parent_cell_id = self.parent_cell_id

        nN = self.config['nNeighbors']
        out = np.nan * np.ones((self.nS, nN))
        for i in range(nN):
            plane_id = cells.plane_id[parent_cell_id[:, i]]
            dens = density.iloc[plane_id].values
            out[:, i] = dens[np.arange(self.nS), self.gene_id]

        return out

