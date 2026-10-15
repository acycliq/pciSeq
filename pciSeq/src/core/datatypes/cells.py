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
        # the cell geometry: a Normal-Inverse-Wishart over the centre and
        # covariance of each cell's gaussian, see the-model/gaussian.md. Starts
        # on the prior: kappa0, nu0, V0 = (nu0 - d - 1) r^2 I, centre at the
        # segmentation centroid. cov keeps E[Sigma] for the ellipses the viewer
        # draws, eig_vals/eig_vecs are of the scale V, which is what the
        # position term in spots.mvn_loglik needs.
        self._d = 3 if config['is3D'] else 2
        kappa0 = config['cell_centroid_prior']['default']
        nu0 = config['cell_cov_prior']['default']
        if nu0 <= self._d + 1:
            raise ValueError('cell_cov_prior must be larger than %d (d + 1) in %dD, got %s'
                             % (self._d + 1, self._d, nu0))
        self._kappa = np.full(self.nC, kappa0, dtype=np.float32)
        self._nu = np.full(self.nC, nu0, dtype=np.float32)
        self._scale = self.ini_scale()
        self._cov = self.ini_cov()
        self._eig_vals, self._eig_vecs = np.linalg.eigh(self._scale)
        self._centroid = self.ini_centroids()
        self._gene_counts = None
        self._gene_counts_prev = None  # the counts before the latest update, see the setter
        self._ini_gene_counts = None  # initial gene counts
        self._background_counts = None
        self._nb_contr = None  # placeholder for the genes' contribution to the negative binomial loglik
        self._theta_bar = None
        self._logtheta_bar = None
        self._nbrs = None
        # the mrf term the last cell_to_cellType used. Diagnostics only, the
        # model reads calc_mrf() rather than this.
        self._mrf = None

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
        # keep the counts we are about to replace. The model never reads this, the
        # diagnostics export does (diagnostics_db.py, scored_gene_count). In the last
        # iteration the class update scores the counts, then spots_to_cell sets them
        # one more time, so at the end of the run these are the counts the class
        # probabilities were worked out from, and geneCount is one update newer.
        # The explanation tools rebuild a cell's score against any class from them
        # (_nb_contributions in mcp/tools.py, nbContributions in the viewer); without
        # them that rebuild is a little off and they refuse to show it. So don't
        # remove it. A reference is enough, geneCount_upd builds a new array each time
        self._gene_counts_prev = self._gene_counts
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
    def mrf(self) -> np.ndarray:
        """The mrf term from the last class update. (nC, nK). Diagnostics only."""
        return self._mrf

    @mrf.setter
    def mrf(self, val):
        self._mrf = val

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
    def dim(self) -> int:
        """2 or 3, the dimensions the cell gaussian lives in."""
        return self._d

    @property
    def cov(self) -> np.ndarray:
        """E[Sigma_c] under the NIW posterior, (nC, d, d). What the ellipses are drawn from."""
        return self._cov

    @cov.setter
    def cov(self, val: np.ndarray):
        self._cov = val

    @property
    def kappa(self) -> np.ndarray:
        """kappa_hat, the pseudo-count behind each cell's centre, (nC,)."""
        return self._kappa

    @kappa.setter
    def kappa(self, val: np.ndarray):
        self._kappa = val

    @property
    def nu(self) -> np.ndarray:
        """nu_hat, the degrees of freedom behind each cell's covariance, (nC,)."""
        return self._nu

    @nu.setter
    def nu(self, val: np.ndarray):
        self._nu = val

    @property
    def scale(self) -> np.ndarray:
        """V_hat, the inverse Wishart scale matrix of each cell, (nC, d, d)."""
        return self._scale

    @scale.setter
    def scale(self, val: np.ndarray):
        self._scale = val

    @property
    def eig_vals(self) -> np.ndarray:
        """Eigenvalues of the scale matrix V_hat, (nC, d)."""
        return self._eig_vals

    @eig_vals.setter
    def eig_vals(self, val: np.ndarray):
        self._eig_vals = val

    @property
    def eig_vecs(self) -> np.ndarray:
        """Eigenvectors of the scale matrix V_hat, (nC, d, d)."""
        return self._eig_vecs

    @eig_vecs.setter
    def eig_vecs(self, val: np.ndarray):
        self._eig_vecs = val

    @property
    def mcr(self) -> float:
        """Returns the mean cell radius."""
        if self.config['cell_radius'] is not None:
            r = self.config['cell_radius']
        else:
            r = self._mcr
        return r

    # per cell, gene and class contributions to the NB loglik from the last class update.
    # check_cell reads it, so dont remove it.
    @property
    def nb_contr(self) -> np.ndarray:
        return self._nb_contr

    @nb_contr.setter
    def nb_contr(self, val):
        self._nb_contr = val

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
        Calculates theta, the per cell scale factor, as a/b.

        Parameters:
            a (np.array): Shape of the Gamma posterior, minus one.
            b (np.array): Rate of the Gamma posterior.
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
        """The prior mean covariance, r^2 I with r the mean cell radius, (nC, d, d)."""
        cov = self.mcr * self.mcr * np.eye(self._d, self._d)
        return np.tile(cov.astype(np.float32), (self.nC, 1, 1))

    def ini_scale(self) -> np.ndarray:
        """V0 = (nu0 - d - 1) r^2 I, the prior scale matrix, (nC, d, d).

        The nu0 - d - 1 is so the prior mean of Sigma is exactly r^2 I, see
        Zhang (2021) on the cell geometry page.
        """
        nu0 = self.config['cell_cov_prior']['default']
        return ((nu0 - self._d - 1) * self.ini_cov()).astype(np.float32)

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

    def nearest_neighbours(self):
        # get the nearest neighbours of each cell
        distances, indices = self.nn().kneighbors(self.zyx_coords)

        # drop the 1st column, it is always the cell itself
        out = {
            'distances': distances[:, 1:],
            'indices': indices[:, 1:]
        }
        return out

    def mrf_support(self):
        """Neighbour support, before mrf_beta.

        Each neighbour's class probs, weighted by how close it is.

        The weight is a gaussian on the distance. The width, sigma, is worked out
        per cell as the median of that cell's own neighbour distances, so what
        counts as close depends on how crowded that bit of tissue is. Note sigma
        is a scale, not a place: you divide the distance by it, you do not measure
        how far the neighbour sits from it.

        Say we are typing cell A and it has two neighbours, B at distance 30 and C
        at distance 60 (the real thing uses nNeighbors of them, usually 9, two is
        just easier to follow)::

            sigma = median(30, 60) = 45
            w_B   = exp(-(30/45)^2) = exp(-0.444) = 0.641
            w_C   = exp(-(60/45)^2) = exp(-1.778) = 0.169

        then rescale so they add up to the number of neighbours::

            w_B -> 1.583,  w_C -> 0.417

        So B's class probabilities go into A's mrf term multiplied by 1.583 and
        C's by 0.417, B counting about 3.8x more than C.

        This is beta-free. calc_mrf below multiplies it by mrf_beta.

        Returns an (nC, nK) array.
        """
        nbrs_idx = self.nbrs['indices']

        # Weight each neighbour by a gaussian on its distance, so closer cells have
        # more influence. The width is per cell, the median of that cell's own
        # neighbour distances, so what counts as close depends on how crowded that
        # bit of tissue is. This is what banksy calls scaled_gaussian.
        #
        # We used to use 1/d here. The trouble with it is that it has no ceiling:
        # when two cells nearly touch, 1/d goes through the roof and that one
        # neighbour walks off with nearly all the weight, up to 7.6 out of 9 on
        # the espio cells. The cell's class then gets decided by a single
        # neighbour and the two of them can sit there flipping each other back and
        # forth forever. exp(-(d/sigma)^2) is flat near zero, so a neighbour that
        # is very close does not gain much by being closer still, and the worst
        # any single one gets is about 3 out of 9.
        #
        # It still weighs by distance, more so than 1/d did if anything: the
        # nearest neighbour counts about 3.2x the furthest, against 2.2x for 1/d,
        # and the weights drop off at every rank for every cell.
        #
        # Normalise so the weights sum to nNeighbors (e.g. 9), matching the
        # scale of zeta (class probs sum to 1 per neighbor, 9 neighbors total).
        # This way proximity and zeta contribute equally to the MRF potential.
        d = self.nbrs['distances']
        nN = nbrs_idx.shape[1]

        sigma = np.median(d, axis=1, keepdims=True)
        nbrs_prxmty = np.exp(-(d / sigma) ** 2)

        # Each row is scaled to sum to nN, so every cell hands out the same total
        # pull no matter how far its neighbours sit.
        nbrs_prxmty = nbrs_prxmty / nbrs_prxmty.sum(axis=1, keepdims=True) * nN

        # Proximity-weighted sum of neighbour class probabilities (zeta)
        nbr_probs = self.classProb[nbrs_idx]  # (nC, nN, nK)
        support = (nbr_probs * nbrs_prxmty[:, :, None]).sum(axis=1)

        # Classes grouped in mrf_pooled_classes share their support, see below
        groups = self.config.get("mrf_pooled_classes") or []
        if groups:
            A = self.class_pooling(groups)
            support = np.einsum('ck, kj -> cj', support, A)

        # If a cell is surrounded by Zero cells, don't let the mrf pull it towards
        # Zero. Zero just never gets any neighbour support. Zero is the last
        # column.
        # 11-Sep-2026: Commenting it out, I have my second thoughts
        # support[:, -1] = 0
        return support

    def class_pooling(self, groups):
        """Makes the classes in a group share their neighbour support.

        Say a cell has 9 neighbours, 8 of them 037 DG Glut and 1 of them
        038 DG-PIR Ex IMN. Normally 037 gets a support of 8 and 038 gets 1, so
        the neighbours push the cell to 037. Put the two in a group and both
        get 8 + 1 = 9. Same bonus for both, so it cancels out and the gene
        counts pick between them.

        This returns a (nK, nK) matrix A and the caller does support @ A. A
        starts as the identity, which leaves the support as it is. Then for
        every group the square block of A for those classes is set to 1, which
        gives each class in the group the sum over the whole group. Classes
        not in any group keep their own support.

        A class can only be in one group. With [class_a, class_b] and
        [class_b, class_c], class_b would get the support of all three while
        class_a and class_c only get two each, so class_b would win just for
        being in both. That raises an error, so does a group with one class
        (it does nothing) and a name that isn't one of the classes.
        """
        class_list = list(self.class_names)
        A = np.eye(len(class_list), dtype=np.float64)
        seen = set()
        for group in groups:
            if len(group) < 2:
                raise ValueError(f"mrf_pooled_classes: {group} needs at least 2 classes")
            missing = [n for n in group if n not in class_list]
            if missing:
                # fail loudly, a typo here would otherwise just do nothing
                raise ValueError(f"mrf_pooled_classes: {missing} not found in the cell classes")
            for n in group:
                if n in seen:
                    raise ValueError(f"mrf_pooled_classes: {n} is in more than one group (or twice in the same one)")
                seen.add(n)
            idx = [class_list.index(n) for n in group]
            A[np.ix_(idx, idx)] = 1
        return A

    def calc_mrf(self):
        """The mrf term: the neighbour support times mrf_beta. (nC, nK).

        Split out from mrf_support because not everything wants a single scalar
        beta. Anything that works out its own per cell or per class coupling
        needs the support on its own, with no beta folded into it yet.
        """
        return self.mrf_support() * self.config["mrf_beta"]

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
