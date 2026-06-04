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

logger = logging.getLogger(__name__)


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
        self._nb_contr = None  # placeholder for the genes' contribution to the negative binomial loglik
        self._mrf = None  # placeholder for the mrf term last used in cell_to_cellType
        self.effective_beta = None  # mrf cap from the last cell_to_cellType call, kept for inspection
        self._theta_bar = None
        self._logtheta_bar = None
        self._nbrs = None

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

    # Property useful only for debugging. Safe to remove
    @property
    def mrf(self) -> np.ndarray:
        return self._mrf

    @mrf.setter
    def mrf(self, val):
        self._mrf = val

    @property
    def ini_gene_counts(self) -> np.ndarray:
        """ Returns an array of shape (nC,) containing the total number of spots
            inside each cell's boundaries.
        """
        return self._ini_gene_counts

    @property
    def theta_bar(self):
        """Returns the eta bar values for genes."""
        return self._theta_bar

    @property
    def logtheta_bar(self):
        """Returns the log eta bar for genes (estimated mean of the posterior)."""
        return self._logtheta_bar

    @property
    def nbrs(self):
        """Returns the nearest neighbors for each cell"""
        return self._nbrs

    @nbrs.setter
    def nbrs(self, val):
        self._nbrs = val

    # -------- METHODS -------- #

    def init_theta(self, a, b):
        """
        Initializes eta values for genes.

        Parameters:
            a (float): Parameter a for eta calculation.
            b (float): Parameter b for eta calculation.
        """
        nK = self.class_names.shape[0]
        self._theta_bar = np.ones([self.nC, nK], dtype=np.float32) * (a / b)
        self._logtheta_bar = np.ones([self.nC, nK], dtype=np.float32) * self._digamma(a, b)

    def calc_theta(self, a, b):
        """
        Calculates eta values for genes.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.
        """
        a = a.astype(np.float32)
        b = b.astype(np.float32)
        self._theta_bar = a[:, None] / b
        # self._logtheta_bar = self._digamma(a, b)

    def _digamma(self, a, b):
        """
        Calculates the digamma function for theta calculation.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.

        Returns:
            np.array: Digamma values.
        """
        return scipy.special.psi(a) - np.log(b)

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

    def nearest_neighbours(self):
        # get the nearest neighbours of each cell
        distances, indices = self.nn().kneighbors(self.zyx_coords)

        # drop the 1st column, it is always the cell itself
        out = {
            'distances': distances[:, 1:],
            'indices': indices[:, 1:]
        }
        return out

    def neighbour_support(self):
        """
        Proximity-weighted soft count of neighbours in each class, after the
        similarity-pooling matrix A. Row c, column k holds n_{c,k}: the
        effective number of neighbours pushing cell c toward class k. This is
        exactly the quantity that beta multiplies in the MRF potential, so
        calc_mrf is just this times beta, and the MRF cap divides by it.

        Returns an (nC, nK) array. The Zero column is identically 0 (the Zero
        row of A is zeroed, so no neighbour ever supports the Zero class).
        """
        nbrs_idx = self.nbrs['indices']

        # Weight each neighbor by 1/distance so closer cells have more influence.
        # Normalise so the weights sum to nNeighbors (e.g. 9), matching the
        # scale of zeta (class probs sum to 1 per neighbor, 9 neighbors total).
        # This way proximity and zeta contribute equally to the MRF potential.
        nbrs_prxmty = 1/self.nbrs['distances']
        nbrs_prxmty = nbrs_prxmty / nbrs_prxmty.sum(axis=1, keepdims=True) * nbrs_idx.shape[1]

        # Proximity-weighted sum of neighbour class probabilities (zeta)
        nbr_probs = self.classProb[nbrs_idx]  # (nC, nN, nK)
        mrf = (nbr_probs * nbrs_prxmty[:, :, None]).sum(axis=1)

        # Row-sum note: at this point sum_k mrf[c, k] = nN per cell. The
        # A-multiplication and the beta scaling below both break this, but
        # for different reasons:
        #   - A modification (Zero-row=0): zeros the Zero column only.
        #     Real-vs-real differences (e.g. Oligo vs Astro) are preserved
        #     exactly, so the softmax over real classes is unchanged.
        #     Zero just loses MRF support.
        #   - beta scaling: multiplies every entry by beta, which SCALES
        #     every class-vs-class difference. This intentionally sharpens
        #     (beta > 1) or flattens (beta < 1) the softmax -- beta is the
        #     parameter that controls how strongly the MRF influences the
        #     cell-class decision.
        # The absolute row sum itself doesn't matter for softmax (which is
        # shift-invariant under adding a constant to every class). What
        # matters is the per-class differences, which A and beta shape on
        # purpose.

        # Similarity matrix A of shape (nK, nK). A[i, j] = 1 means a neighbour
        # classified as class j contributes to the MRF support of class i (rows
        # are receivers, columns are donors). The identity diagonal is the
        # standard case: a neighbour of class k supports the cell under focus
        # being class k, and contributes nothing to any other class.
        # Symmetric off-diagonal 1s pool two similar classes: a neighbour of
        # either class supports both, which neutralises the MRF between them and
        # leaves the gene log-likelihood to pick the winner. Without this, a
        # rare class (e.g. 038 DG-PIR Ex IMN) embedded inside a dense majority
        # (037 DG Glut) loses the softmax to its sister class even when the gene
        # evidence slightly favours it, because the neighbourhood votes are
        # overwhelmingly for the majority. See notes/mrf_similarity_matrix.md
        # for the full derivation and a worked toy example.
        class_list = list(self.class_names)
        A = np.eye(len(class_list), dtype=mrf.dtype)
        for a, b in self.config["similarity_pairs"]:
            if a in class_list and b in class_list:
                ia, ib = class_list.index(a), class_list.index(b)
                A[ia, ib] = A[ib, ia] = 1

        # Zero-classified neighbours contribute no MRF support to any class.
        # Without this, a cell surrounded by Zero neighbours gets dragged toward
        # Zero by neighbour pressure, which we don't want -- the data should
        # decide whether the cell is Zero, not the neighbourhood. This breaks
        # the symmetry of A (Zero column is unchanged, Zero row is now all zeros).
        assert class_list[-1] == 'Zero', "Last class must be Zero"
        A[-1, :] = 0

        mrf = oe.contract('ck, kj -> cj', mrf, A)

        return mrf

    def calc_mrf(self, effective_beta=None):
        mrf = self.neighbour_support()

        if effective_beta is None:
            out = mrf * self.config["mrf_beta"]
        else:
            out = mrf * effective_beta

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
