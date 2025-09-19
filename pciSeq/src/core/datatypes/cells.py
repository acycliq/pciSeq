# Standard library imports
import logging
from typing import Tuple, Dict, Any

# Third party imports
import numpy as np
import pandas as pd
import scipy
from natsort import natsort_keygen
from sklearn.neighbors import NearestNeighbors
import numpy_groupies as npg
import opt_einsum as oe

# Local imports
from ..utils.cell_utils import read_image_objects, keep_labels_unique
from ..utils import ops_utils as utils
from ..utils.geometry import anisotropy_calc

cells_logger = logging.getLogger(__name__)


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
        self._eig_vals, self._eig_vecs = np.linalg.eigh(self._cov)
        self.nu_0 = config['mean_gene_counts_per_cell']
        self._centroid = self.ini_centroids()
        self._gene_counts = None
        self._ini_gene_counts = None  # initial gene counts
        self._background_counts = None
        self.on_planes = dict(zip(_cells_df['label'], _cells_df['values']))
        self._nb_contr = None  # placeholder for the genes' contribution to the negative binomial loglik
        self._plane_adj = None

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
    def eig_vals(self) -> np.ndarray:
        """Returns the eigenvalues of the covariance matrix."""
        return self._eig_vals

    @eig_vals.setter
    def eig_vals(self, val: np.ndarray):
        """Sets the eigenvalues of the covariance matrix."""
        self._eig_vals = val

    @property
    def eig_vecs(self) -> np.ndarray:
        """Returns the eigenvectors of the covariance matrix."""
        return self._eig_vecs

    @eig_vecs.setter
    def eig_vecs(self, val: np.ndarray):
        """Sets the eigenvectors of the covariance matrix."""
        self._eig_vecs = val

    @property
    def mcr(self) -> float:
        """Returns the mean cell radius."""
        if self.config['cell_radius'] is not None:
            r = self.config['cell_radius']
        else:
            r = self._mcr
        return r

    # Property useful only for debugging. Safe to remove
    @property
    def nb_contr(self) -> np.ndarray:
        return self._nb_contr

    @nb_contr.setter
    def nb_contr(self, val):
        self._nb_contr = val

    @property
    def plane_id(self) -> np.ndarray:
        cell_coords = anisotropy_calc(self.centroid.values, voxel_size=self.config['voxel_size'], inverse=True)
        return np.floor(cell_coords[:,-1]).astype(np.int32)

    @property
    def plane_adj(self) -> pd.DataFrame:
        return self._plane_adj


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
        dim = 3 if self.config['is3D'] else 2
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
        _id = spots.parent_cell_id[:, :-1]
        xyz_spots = spots.xyz_coords
        # out = self.ini_cov() * self.nu_0
        out = np.zeros(self.ini_cov().shape)

        mu_x = mu_bar[_id, 0]  # array of size [nS, N] with the x-coord of the centroid of the N closest cells
        mu_y = mu_bar[_id, 1]  # array of size [nS, N] with the y-coord of the centroid of the N closest cells
        mu_z = mu_bar[_id, 2]  # array of size [nS, N] with the z-coord of the centroid of the N closest cells

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
        agg_00 = npg.aggregate(_id.ravel(), el_00.ravel(), size=self.nC)
        agg_11 = npg.aggregate(_id.ravel(), el_11.ravel(), size=self.nC)
        agg_22 = npg.aggregate(_id.ravel(), el_22.ravel(), size=self.nC)

        agg_01 = npg.aggregate(_id.ravel(), el_01.ravel(), size=self.nC)
        agg_02 = npg.aggregate(_id.ravel(), el_02.ravel(), size=self.nC)
        agg_12 = npg.aggregate(_id.ravel(), el_12.ravel(), size=self.nC)

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


    def calc_plane_adj(self, spots, config):
        """
        Returns an array of shape (n_planes, n_genes) where each row
        represents a plane and each column represents a gene. The first row
        corresponds to the background plane (with negative plane_id) and contains
        all 1.0 values. Subsequent rows correspond to image planes and contain
        computed adjustment factors for the single cell data.
        """

        background_plane_id = self.plane_id[0]  # negative index for background
        regular_plane_ids = self.plane_id[1:]   # actual plane indices

        density = utils.gene_density(spots, config)
        out = density.iloc[regular_plane_ids]

        # Add background row (all 1s - no adjustment)
        background_row = pd.DataFrame(1.0, index=[background_plane_id], columns=out.columns)
        out = pd.concat([background_row, out], axis=0)

        return out


    # -------------------------- CONVENIENCE METHODS ----------------------- #
    def gene_reads_per_class(self):
        """Calculate total (weighted by class prob) gene reads for each class.

        Returns:
            np.ndarray: Shape (G, K) total reads per class and gene
        """
        # Calculate weighted sum of gene reads for each class and gene using classProb as weights
        weighted_sum = oe.contract('cg, ck -> gk', self.geneCount, self.classProb, optimize='optimal')
        return weighted_sum

    def mean_gene_reads_per_class(self):
        """Calculate the average gene reads for each cell class/type in a soft clustering setup.

        In soft clustering, each cell belongs to multiple classes with probabilities \( w_{ck} \).
        The average number of reads for gene \( g \) in class \( k \) is computed as:

        \[
        \overline{r}_{gk} = \frac{\sum_{c=1}^{C} x_{cg} \cdot w_{ck}}{\sum_{c=1}^{C} w_{ck}}
        \]

        Where:
            - \( x_{cg} \): Number of reads for gene \( g \) in cell \( c \)
            - \( w_{ck} \): Probability that cell \( c \) belongs to class \( k \)
            - The numerator is the total weighted sum of reads for gene \( g \) in class \( k \)
            - The denominator is the total probability mass of class \( k \)

        Returns:
            np.ndarray: Shape (G, K), where:
                G = number of genes
                K = number of cell classes/types
        """
        weighted_sum = self.gene_reads_per_class()

        # Calculate total probability mass (size) for each class
        class_totals = self.classProb.sum(axis=0)

        # Assuming weighted_sum has shape (319, 39) and class_totals has shape (39,)
        result = np.divide(
            weighted_sum,
            class_totals,
            out=np.zeros_like(weighted_sum),
            where=class_totals != 0
        )

        return result
