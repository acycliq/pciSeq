"""
Core Data Types Module for pciSeq

This module defines the fundamental data structures used throughout pciSeq, a package designed for
analyzing and processing image segmentation data, particularly in the context of cell and gene analysis.

Key Functionalities:
- Manages and processes cell segmentation data, gene expression data, and RNA spot detection.
- Provides structures for handling single-cell RNA sequencing data and classifying cell types.

Key Classes:
-----------
Cells:
    Handles cell segmentation data, including properties like centroids,
    covariance matrices, and gene counts. It provides methods for calculating
    nearest neighbors and scatter matrices, essential for spatial analysis.

Genes:
    Manages gene-specific data and calculations, including initialization and
    computation of gene expression parameters. This class is crucial for
    understanding gene distribution and expression levels across cells.

Spots:
    Processes RNA spot detection data, including spatial coordinates and
    cell assignments. It includes methods for reading spot data, calculating
    nearest cells, and managing gene counts, which are vital for spatial
    transcriptomics.

SingleCell:
    Handles single-cell RNA sequencing reference data, including mean
    expression levels per cell type. This class supports the integration of
    single-cell data into broader analyses, providing a reference for gene
    expression.

CellType:
    Manages cell type classification, including prior probabilities and
    class assignments. It helps in understanding the distribution of different
    cell types and their characteristics within a dataset.

Notes:
------
- All numerical computations use numpy arrays for efficiency, ensuring fast
  processing of large datasets.
- Sparse matrices are used where appropriate for memory efficiency, particularly
  in handling large, sparse datasets.
- Dask is used for delayed computations of large arrays, allowing for scalable
  and parallel processing.

Dependencies:
------------
- scipy: For statistical computations and special functions.
- numpy: For efficient numerical operations and array handling.
- pandas: For data management and manipulation, particularly with tabular data.
- dask: For delayed computations and handling large datasets efficiently.
- sklearn: For nearest neighbor calculations, essential for spatial analysis.
- numpy_groupies: For efficient grouping operations, used in aggregating data.

This module is designed to be integrated into larger pipelines for image segmentation
and analysis, providing robust data structures and methods for handling complex
biological data.
"""

import scipy
import numpy as np
import pandas as pd
import numpy_groupies as npg
from natsort import natsort_keygen
from .utils import read_image_objects, keep_labels_unique
from sklearn.neighbors import NearestNeighbors
import logging
from typing import Tuple, Dict, Any

datatypes_logger = logging.getLogger(__name__)


class Cells(object):
    """
    Represents cell segmentation data, including properties like centroids,
    covariance matrices, and gene counts. Provides methods for calculating
    nearest neighbors and scatter matrices.

    Attributes:
        config (dict): Configuration parameters for cell data.
        ini_cell_props (dict): Initial cell properties.
        nC (int): Number of cells.
        classProb (np.array): Class probabilities for cells.
        class_names (list): Names of cell classes.
        _cov (np.array): Covariance matrices for cells.
        nu_0 (float): Mean gene counts per cell.
        _centroid (pd.DataFrame): Centroid coordinates for cells.
        _gene_counts (np.array): Gene counts for cells.
        _background_counts (np.array): Num of spots assigned to background.
    """

    def __init__(self, _cells_df: pd.DataFrame, config: Dict):
        """
        Initializes the Cells object with cell data and configuration.

        Parameters:
            _cells_df (pd.DataFrame): DataFrame containing cell data.
            config (dict): Configuration parameters for cell data.
        """
        self.config = config
        self.ini_cell_props, self._mcr = read_image_objects(_cells_df, config)
        self.nC = len(self.ini_cell_props['cell_label'])
        self.classProb = None
        self.class_names = None
        self._cov = self.ini_cov()
        self.nu_0 = config['mean_gene_counts_per_cell']
        self._centroid = self.ini_centroids()
        self._gene_counts = None
        self._background_counts = None

    # -------- PROPERTIES -------- #
    @property
    def zyx_coords(self) -> np.ndarray:
        """Returns the centroid coordinates in z, y, x order."""
        # Convenience property but maybe it should be removed. Potentially could be using memory
        # lots of memory for no real benefit

        return self.centroid[['z', 'y', 'x']].values

    @property
    def geneCount(self) -> np.ndarray:
        """Returns the gene counts for cells."""
        return self._gene_counts

    @geneCount.setter
    def geneCount(self, val: np.ndarray):
        """Sets the gene counts for cells."""
        self._gene_counts = val

    @property
    def background_counts(self) -> np.ndarray:
        """Returns the background counts for cells."""
        return self._background_counts

    @background_counts.setter
    def background_counts(self, val):
        """Sets the background gene counts."""
        self._background_counts = val

    @property
    def total_counts(self) -> np.ndarray:
        """Returns the total gene counts for cells."""
        return self.geneCount.sum(axis=1)

    @property
    def centroid(self) -> pd.DataFrame:
        """Returns a copy of the centroid DataFrame."""
        return self._centroid.copy()

    @centroid.setter
    def centroid(self, df: pd.DataFrame):
        """
        Sets the centroid DataFrame.

        Parameters:
            df (pd.DataFrame): DataFrame containing centroid coordinates.
        """
        assert isinstance(df, pd.DataFrame), 'Input should be a dataframe'
        assert set(df.columns.values) == {'x', 'y', 'z'}, 'Dataframe columns should be ''x'', ''y'' and ''z'' '
        df.index.name = 'cell_label'
        self._centroid = df.copy()

    @property
    def cov(self) -> np.ndarray:
        """Returns the covariance matrices for cells."""
        return self._cov

    @cov.setter
    def cov(self, val: np.ndarray):
        """Sets the covariance matrices for cells."""
        self._cov = val

    @property
    def mcr(self) -> float:
        """Returns the mean cell radius."""
        if self.config['cell_radius'] is not None:
            r = self.config['cell_radius']
        else:
            r = self._mcr
        return r

    # -------- METHODS -------- #
    def ini_centroids(self) -> pd.DataFrame:
        """
        Initializes the centroids for cells.

        Returns:
            pd.DataFrame: DataFrame containing centroid coordinates.
        """
        d = {
            'x': self.ini_cell_props['x0'],
            'y': self.ini_cell_props['y0'],
            'z': self.ini_cell_props['z0'],
        }
        df = pd.DataFrame(d)
        return df.copy()

    def ini_cov(self) -> np.ndarray:
        """
        Initializes the covariance matrices for cells.

        Returns:
            np.array: Array of covariance matrices.
        """
        dim = 3
        cov = self.mcr * self.mcr * np.eye(dim, dim)
        return np.tile(cov.astype(np.float32), (self.nC, 1, 1))

    def nn(self) -> NearestNeighbors:
        """
        Calculates the nearest neighbors for cells.

        Returns:
            NearestNeighbors: Fitted NearestNeighbors object.
        """
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.zyx_coords)
        return nbrs

    def scatter_matrix(self, spots: 'Spots') -> np.ndarray:
        """
        Calculates the scatter matrix for cells based on spot data.

        Parameters:
            spots (Spots): Spots object containing spot data.

        Returns:
            np.array: Scatter matrix for cells.
        """
        mu_bar = self.centroid.values
        prob = spots.parent_cell_prob[:, :-1]
        id = spots.parent_cell_id[:, :-1]
        xyz_spots = spots.xyz_coords
        out = self.ini_cov() * self.nu_0

        mu_x = mu_bar[id, 0]  # array of size [nS, N] with the x-coord of the centroid of the N closest cells
        mu_y = mu_bar[id, 1]  # array of size [nS, N] with the y-coord of the centroid of the N closest cells
        mu_z = mu_bar[id, 2]  # array of size [nS, N] with the z-coord of the centroid of the N closest cells

        N = mu_x.shape[1]
        _x = np.tile(xyz_spots[:, 0], (N, 1)).T  # array of size [nS, N] populated with the x-coord of the spot
        _y = np.tile(xyz_spots[:, 1], (N, 1)).T  # array of size [nS, N] populated with the y-coord of the spot
        _z = np.tile(xyz_spots[:, 2], (N, 1)).T  # array of size [nS, N] populated with the z-coord of the spot

        x_centered = _x - mu_x  # subtract the cell centroid x-coord from the spot x-coord
        y_centered = _y - mu_y  # subtract the cell centroid y-coord from the spot y-coord
        z_centered = _z - mu_z  # subtract the cell centroid z-coord from the spot z-coord

        el_00 = prob * x_centered * x_centered  # contribution to the scatter matrix's [0, 0] element
        el_11 = prob * y_centered * y_centered  # contribution to the scatter matrix's [1, 1] element
        el_22 = prob * z_centered * z_centered  # contribution to the scatter matrix's [2, 2] element

        el_01 = prob * x_centered * y_centered  # contribution to the scatter matrix's [0, 1] element
        el_02 = prob * x_centered * z_centered  # contribution to the scatter matrix's [0, 2] element
        el_12 = prob * y_centered * z_centered  # contribution to the scatter matrix's [1, 2] element

        # Aggregate all contributions to get the scatter matrix
        agg_00 = npg.aggregate(id.ravel(), el_00.ravel(), size=self.nC)
        agg_11 = npg.aggregate(id.ravel(), el_11.ravel(), size=self.nC)
        agg_22 = npg.aggregate(id.ravel(), el_22.ravel(), size=self.nC)

        agg_01 = npg.aggregate(id.ravel(), el_01.ravel(), size=self.nC)
        agg_02 = npg.aggregate(id.ravel(), el_02.ravel(), size=self.nC)
        agg_12 = npg.aggregate(id.ravel(), el_12.ravel(), size=self.nC)

        # Return now the scatter matrix. Some cell might not have any spots nearby. For those empty cells,
        # the scatter matrix will be a squared zero array. That is fine.
        out[:, 0, 0] = agg_00
        out[:, 1, 1] = agg_11
        out[:, 2, 2] = agg_22

        out[:, 0, 1] = agg_01
        out[:, 0, 2] = agg_02
        out[:, 1, 2] = agg_12

        out[:, 1, 0] = agg_01
        out[:, 2, 0] = agg_02
        out[:, 2, 1] = agg_12

        return out.astype(np.float32)


# ----------------------------------------Class: Genes--------------------------------------------------- #
class Genes(object):
    """
    Manages gene-specific data and calculations, including initialization and
    computation of gene expression parameters.

    Attributes:
        gene_panel (np.array): Array of unique gene names.
        _eta_bar (np.array): Eta bar values: This is basically the expected Gene inefficiency.
        _logeta_bar (np.array): Log eta bar values for genes.
        nG (int): Number of genes.
    """

    def __init__(self, spots):
        """
        Initializes the Genes object with spot data.

        Parameters:
            spots (Spots): Spots object containing spot data.
        """
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta_bar = None
        self._logeta_bar = None
        self.nG = len(self.gene_panel)

    @property
    def eta_bar(self):
        """Returns the eta bar values for genes."""
        return self._eta_bar

    @property
    def logeta_bar(self):
        """Returns the log eta bar values for genes."""
        return self._logeta_bar

    def init_eta(self, a, b):
        """
        Initializes eta values for genes.

        Parameters:
            a (float): Parameter a for eta calculation.
            b (float): Parameter b for eta calculation.
        """
        self._eta_bar = np.ones(self.nG, dtype=np.float32) * (a / b)
        self._logeta_bar = np.ones(self.nG, dtype=np.float32) * self._digamma(a, b)

    def calc_eta(self, a, b):
        """
        Calculates eta values for genes.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.
        """
        a = a.astype(np.float32)
        b = b.astype(np.float32)
        self._eta_bar = a / b
        self._logeta_bar = self._digamma(a, b)

    def _digamma(self, a, b):
        """
        Calculates the digamma function for eta calculation.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.

        Returns:
            np.array: Digamma values.
        """
        return scipy.special.psi(a) - np.log(b)


# ----------------------------------------Class: Spots--------------------------------------------------- #
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

    # -------- PROPERTIES -------- #
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
        lst = list(zip(*[self.data.x, self.data.y, self.data.z]))
        return np.array(lst, dtype=np.float32)

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

    # -------- METHODS -------- #
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
        if ~is3D:
            # that shouldn't really be necessary. If the data are 2d then the z dimension if just a dummy dimension.
            # and inference should still hold (with the dummy z dimension)
            # I am just removing z here for backwards compatibility; To yield the same results as in the 2d case.
            # Otherwise, the dummy z dimension will get intp the loglikelihood calculations and the results will
            # differ. Still though they should be regarded as valid results
            data = data[:, :-1]
            centroids = centroids[:, :-1]
            covs = covs[:, :-1, :-1]
        out = self.multiple_logpdfs(data, centroids, covs)
        return out

    def multiple_logpdfs(self, x: np.ndarray, means: np.ndarray, covs: np.ndarray) -> np.ndarray:
        """
        vectorised mvn log likelihood evaluated at multiple pairs of (centroid_1, cov_1), ..., (centroid_N, cov_N)
        Taken from http://gregorygundersen.com/blog/2020/12/12/group-multivariate-normal-pdf/
        """
        # Thankfully, NumPy broadcasts `eigh`.
        vals, vecs = np.linalg.eigh(covs)

        # Compute the log determinants across the second axis.
        logdets = np.sum(np.log(vals), axis=1)

        # Invert the eigenvalues.
        valsinvs = 1. / vals

        # Add a dimension to `valsinvs` so that NumPy broadcasts appropriately.
        Us = vecs * np.sqrt(valsinvs)[:, None]
        devs = x - means

        # Use `einsum` for matrix-vector multiplications across the first dimension.
        devUs = np.einsum('ni,nij->nj', devs, Us)

        # Compute the Mahalanobis distance by squaring each term and summing.
        mahas = np.sum(np.square(devUs), axis=1)

        # Compute and broadcast scalar normalizers.
        dim = len(vals[0])
        log2pi = np.log(2 * np.pi)

        return -0.5 * (dim * log2pi + mahas + logdets)

    def zero_class_counts(self, geneNo, pCellZero):
        """
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
        Calculates the expectation of a gamma distribution.

        Parameters:
            rho (np.array): Shape parameters.
            beta (np.array): Rate parameters.

        Returns:
            np.array: Expected values.
        """
        r = rho[:, :, None]
        return r / beta

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


# ----------------------------------------Class: SingleCell--------------------------------------------------- #
class SingleCell(object):
    """
    Handles single-cell RNA sequencing reference data, including mean
    expression levels per cell type. Supports integration of single-cell
    data into broader analyses.

    Attributes:
        isMissing (bool): Indicates if single-cell data is missing.
        config (dict): Configuration parameters for single-cell data.
        _mean_expression (pd.DataFrame): Mean expression levels.
        _log_mean_expression (pd.DataFrame): Log mean expression levels.
    """

    def __init__(self, scdata: pd.DataFrame, genes: np.ndarray, config: Dict):
        """
        Initializes the SingleCell object with single-cell data and configuration.

        Parameters:
            scdata (pd.DataFrame): Single-cell data.
            genes (np.array): Array of gene names.
            config (dict): Configuration parameters for single-cell data.
        """
        self.isMissing = None  # Will be set to False if single cell data are assumed known and given as an input
        # otherwise, if they are unknown, this will be set to True and the algorithm will
        # try to estimate them
        # self.raw_data = self._raw_data(scdata, genes)
        self.config = config
        self._mean_expression, self._log_mean_expression = self._setup(scdata, genes, self.config)

    def _setup(self, scdata: pd.DataFrame, genes: np.ndarray, config: Dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Sets up the mean and log mean expression levels.

        Parameters:
            scdata (pd.DataFrame): Single-cell data.
            genes (np.array): Array of gene names.
            config (dict): Configuration parameters.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Mean and log mean expression levels.
        """
        if scdata is None:
            datatypes_logger.info('Single Cell data are missing. Cannot determine mean expression per cell class.')
            datatypes_logger.info('We will try to estimate the array instead')
            datatypes_logger.info('Starting point is a diagonal array of size numGenes-by-numGenes')
            expr = self._diag(genes)
            self.isMissing = True
        else:
            expr = self._raw_data(scdata, genes)
            self.isMissing = False

        self.raw_data = expr
        me, lme = self._helper(expr.copy())

        assert me.columns[-1] == 'Zero', "Last column should be the Zero class"
        assert lme.columns[-1] == 'Zero', "Last column should be the Zero class"
        return me.astype(np.float32), lme.astype(np.float32)

    # -------- PROPERTIES -------- #
    @property
    def mean_expression(self):
        """Returns the mean expression levels."""
        assert self._mean_expression.columns[-1] == 'Zero', "Last column should be the Zero class"
        return self._mean_expression

    @property
    def log_mean_expression(self):
        """Returns the log mean expression levels."""
        assert self._log_mean_expression.columns[-1] == 'Zero', "Last column should be the Zero class"
        return self._log_mean_expression

    @property
    def genes(self):
        """Returns the gene names."""
        return self.mean_expression.index.values

    @property
    def classes(self):
        """Returns the class names."""
        return self.mean_expression.columns.values

    ## Helper functions ##
    def _set_axes(self, df):
        """
        Sets the axes labels for a DataFrame.

        Parameters:
            df (pd.DataFrame): DataFrame to set axes for.

        Returns:
            pd.DataFrame: DataFrame with set axes.
        """
        df = df.rename_axis("class_name", axis="columns").rename_axis('gene_name')
        return df

    def _remove_zero_cols(self, df):
        """
        Removes columns with all zero values from a DataFrame.

        Parameters:
            df (pd.DataFrame): DataFrame to remove zero columns from.

        Returns:
            pd.DataFrame: DataFrame with zero columns removed.
        """
        out = df.loc[:, (df != 0).any(axis=0)]
        return out

    def _helper(self, expr):
        """
        Helper function to process expression data.

        Parameters:
            expr (pd.DataFrame): Expression data.

        Returns:
            tuple: Processed mean and log-mean expression data.
        """

        # order by column name
        expr = expr.copy().sort_index(axis=0).sort_index(axis=1, key=natsort_keygen(key=lambda y: y.str.lower()))

        # append at the end the Zero class
        expr['Zero'] = np.zeros([expr.shape[0], 1])
        me = expr.rename_axis('gene_name').rename_axis("class_name", axis="columns")

        # add the regularization parameter
        me = me + self.config['SpotReg']

        # log mean expression
        lme = np.log(me)
        return me, lme

    def _gene_expressions(self, fitted, scale):
        """
        Calculates expected mean gene counts. The prior *IS NOT* taken
        into account. We use data evidence only
        For the zero class only the prior is used *AND NOT* data
        evidence.

        Parameters:
            fitted (np.array): Fitted values.
            scale (np.array): Scale values.

        Returns:
            tuple: Mean and log-mean gene expressions.
        """

        # the prior on mean expression follows a Gamma(m * M , m), where M is the starting point (the initial
        # array) of single cell data
        # 07-May-2023. Hiding m from the config.py. Should bring it back at a later version
        # m = self.config['m']
        m = 1
        a = fitted + m * (self.raw_data + self.config['SpotReg'])
        b = scale + m
        me = a / b
        lme = scipy.special.psi(a) - np.log(b)

        # the expressions for the zero class are a 0.0 plus the regularition param
        zero_col = np.zeros(me.shape[0]) + self.config['SpotReg']
        me = me.assign(Zero=zero_col)

        # For the mean of the log-expressions, again only the prior is used for the Zero class
        zero_col_2 = scipy.special.psi(m * zero_col) - np.log(m)
        lme = lme.assign(Zero=zero_col_2)
        return me, lme

    def _raw_data(self, scdata: pd.DataFrame, genes: np.ndarray) -> pd.DataFrame:
        """
        Processes raw single-cell data, filtering out any genes outside the gene panel and grouping by cell type.

        Parameters:
            scdata (pd.DataFrame): Single-cell data.
            genes (np.array): Array of gene names.

        Returns:
            pd.DataFrame: Processed single-cell data.
        """
        assert np.all(scdata >= 0), "Single cell dataframe has negative values"
        datatypes_logger.info('Single cell data passed-in have %d genes and %d cells' % (scdata.shape[0], scdata.shape[1]))
        datatypes_logger.info('Single cell data: Keeping counts for the gene panel of %d only' % len(genes))
        df = scdata.loc[genes]

        # set the axes labels
        df = self._set_axes(df)

        # remove any rows with the same gene label
        df = keep_labels_unique(df)

        df = self._remove_zero_cols(df.copy())
        dfT = df.T

        datatypes_logger.info('Single cell data: Grouping gene counts by cell type. Aggregating function is the mean.')
        out = dfT.groupby(dfT.index.values).agg('mean').T
        datatypes_logger.info('Grouped single cell data have %d genes and %d cell types' % (out.shape[0], out.shape[1]))
        return out

    def _diag(self, genes):
        """
        Creates a diagonal matrix for single-cell data initialization.

        Parameters:
            genes (np.array): Array of gene names.

        Returns:
            pd.DataFrame: Diagonal single-cell data.
        """
        nG = len(genes)
        mgc = self.config['mean_gene_counts_per_class']
        arr = mgc * np.eye(nG)
        labels = ['class_%d' % (i + 1) for i, _ in enumerate(genes)]
        df = pd.DataFrame(arr).set_index(genes)
        df.columns = labels
        return df


# ---------------------------------------- Class: CellType --------------------------------------------------- #
class CellType(object):
    """
    Manages cell type classification, including prior probabilities and
    class assignments. Helps in understanding the distribution of different
    cell types and their characteristics within a dataset.

    Attributes:
        _names (np.array): Names of cell types.
        _alpha (np.array): Alpha values for cell types.
        config (dict): Configuration parameters for cell types.
        single_cell_data_missing (bool): Indicates if single-cell data is missing.
    """

    def __init__(self, single_cell: SingleCell, config: Dict):
        """
        Initializes the CellType object with single-cell data and configuration.

        Parameters:
            single_cell (SingleCell): SingleCell object containing single-cell data.
            config (dict): Configuration parameters for cell types.
        """
        assert single_cell.classes[-1] == 'Zero', "Last label should be the Zero class"
        self._names = single_cell.classes
        self._alpha = None
        self.config = config
        self.single_cell_data_missing = single_cell.isMissing

    @property
    def names(self) -> np.ndarray:
        """Returns the names of cell types."""
        assert self._names[-1] == 'Zero', "Last label should be the Zero class"
        return self._names

    @property
    def nK(self) -> int:
        """Returns the number of cell types."""
        return len(self.names)

    @property
    def alpha(self) -> np.ndarray:
        """Returns the alpha values for cell types."""
        return self._alpha

    @alpha.setter
    def alpha(self, val: np.ndarray):
        """Sets the alpha values for cell types."""
        self._alpha = val

    @property
    def pi_bar(self) -> np.ndarray:
        """Returns the pi bar values for cell types."""
        return self.alpha / self.alpha.sum()

    @property
    def logpi_bar(self) -> np.ndarray:
        """Returns the log pi bar values for cell types."""
        return scipy.special.psi(self.alpha) - scipy.special.psi(self.alpha.sum())

    @property
    def prior(self) -> np.ndarray:
        """Returns the prior probabilities for cell types."""
        return self.pi_bar

    @property
    def log_prior(self) -> np.ndarray:
        """Returns the log prior probabilities for cell types."""
        if self.single_cell_data_missing or self.config['cell_type_prior'] == 'weighted':
            return self.logpi_bar
        else:
            return np.log(self.prior)

    def size(self, cells: Cells) -> np.ndarray:
        """
        Calculates the size of each cell type, i.e., the number of cells in each type.

        Parameters:
            cells (Cells): Cells object containing cell data.

        Returns:
            np.array: Sizes of cell types.
        """
        return cells.classProb.sum(axis=0)

    def ini_prior(self):
        """Initializes the prior probabilities for cell types."""
        self.alpha = self.ini_alpha()

    def ini_alpha(self) -> np.ndarray:
        """
        Initializes the alpha values for cell types.

        Returns:
            np.array: Initialized alpha values.
        """
        ones = np.ones(self.nK - 1)
        return np.append(ones, sum(ones)).astype(np.float32)
