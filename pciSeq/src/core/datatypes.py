import scipy
import numpy as np
import pandas as pd
import dask
import numpy_groupies as npg
from natsort import natsort_keygen
from .utils import read_image_objects, keep_labels_unique
from sklearn.neighbors import NearestNeighbors
import logging

datatypes_logger = logging.getLogger(__name__)


class Cells(object):
    # Get rid of the properties where not necessary!!
    def __init__(self, _cells_df, config):
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
    def yx_coords(self):
        return self.centroid[['y', 'x']].values

    @property
    def geneCount(self):
        return self._gene_counts

    @geneCount.setter
    def geneCount(self, val):
        self._gene_counts = val

    @property
    def background_counts(self):
        return self._background_counts

    @background_counts.setter
    def background_counts(self, val):
        # assert val[1:, :].sum() == 0, 'Input array must be zero everywhere apart from the top row'
        # self._background_counts = val[0, :]
        self._background_counts = val

    @property
    def total_counts(self):
        # tc = self.geneCount.sum(axis=1)
        return self.geneCount.sum(axis=1)

    @property
    def centroid(self):
        # lst = list(zip(*[self._centroid['x'], self._centroid['y']]))
        return self._centroid.copy()

    @centroid.setter
    def centroid(self, df):
        assert isinstance(df, pd.DataFrame), 'Input should be a dataframe'
        assert set(df.columns.values) == {'x', 'y'}, 'Dataframe columns should be ''x'', ''y'' '
        df.index.name = 'cell_label'
        self._centroid = df.copy()

    @property
    def cov(self):
        return self._cov

    @cov.setter
    def cov(self, val):
        self._cov = val

    @property
    def mcr(self):
        if self.config['cell_radius'] is not None:
            r = self.config['cell_radius']
        else:
            r = self._mcr
        return r

    # -------- METHODS -------- #
    def ini_centroids(self):
        d = {
            'x': self.ini_cell_props['x0'],
            'y': self.ini_cell_props['y0'],
        }
        df = pd.DataFrame(d)
        return df.copy()

    def ini_cov(self):
        mcr = self.mcr
        cov = mcr * mcr * np.eye(2, 2, dtype=np.float32)
        return np.tile(cov, (self.nC, 1, 1))

    # def dapi_mean_cell_radius(self):
    #     return np.nanmean(np.sqrt(self.ini_cell_props['area'] / np.pi)) * 0.5

    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.yx_coords)
        return nbrs

    def scatter_matrix(self, spots):
        mu_bar = self.centroid.values
        prob = spots.parent_cell_prob[:, :-1]
        id = spots.parent_cell_id[:, :-1]
        xy_spots = spots.xy_coords
        out = self.ini_cov() * self.nu_0

        mu_x = mu_bar[id, 0]  # array of size [nS, N] with the x-coord of the centroid of the N closest cells
        mu_y = mu_bar[id, 1]  # array of size [nS, N] with the y-coord of the centroid of the N closest cells

        N = mu_x.shape[1]
        _x = np.tile(xy_spots[:, 0], (N, 1)).T  # array of size [nS, N] populated with the x-coord of the spot
        _y = np.tile(xy_spots[:, 1], (N, 1)).T  # array of size [nS, N] populated with the y-coord of the spot
        x_centered = _x - mu_x  # subtract the cell centroid x-coord from the spot x-coord
        y_centered = _y - mu_y  # subtract the cell centroid y-coord from the spot x-coord

        el_00 = prob * x_centered * x_centered  # contribution to the scatter matrix's [0, 0] element
        el_11 = prob * y_centered * y_centered  # contribution to the scatter matrix's off-diagonal element
        el_01 = prob * x_centered * y_centered  # contribution to the scatter matrix's [1, 1] element

        # Aggregate all contributions to get the scatter matrix
        agg_00 = npg.aggregate(id.ravel(), el_00.ravel(), size=self.nC)
        agg_11 = npg.aggregate(id.ravel(), el_11.ravel(), size=self.nC)
        agg_01 = npg.aggregate(id.ravel(), el_01.ravel(), size=self.nC)

        # Return now the scatter matrix. Some cell might not have any spots nearby. For those empty cells,
        # the scatter matrix will be a squared zero array. That is fine.
        out[:, 0, 0] = agg_00
        out[:, 1, 1] = agg_11
        out[:, 0, 1] = agg_01
        out[:, 1, 0] = agg_01

        return out.astype(np.float32)


# ----------------------------------------Class: Genes--------------------------------------------------- #
class Genes(object):
    def __init__(self, spots):
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta_bar = None
        self._logeta_bar = None
        self.nG = len(self.gene_panel)

    @property
    def eta_bar(self):
        return self._eta_bar

    @property
    def logeta_bar(self):
        return self._logeta_bar

    def init_eta(self, a, b):
        self._eta_bar = np.ones(self.nG, dtype=np.float32) * (a / b)
        self._logeta_bar = np.ones(self.nG, dtype=np.float32) * self._digamma(a, b)

    def calc_eta(self, a, b):
        a = a.astype(np.float32)
        b = b.astype(np.float32)
        self._eta_bar = a / b
        self._logeta_bar = self._digamma(a, b)

    def _digamma(self, a, b):
        return scipy.special.psi(a) - np.log(b)


# ----------------------------------------Class: Spots--------------------------------------------------- #
class Spots(object):
    def __init__(self, spots_df, config):
        self._parent_cell_prob = None
        self._parent_cell_id = None
        self.Dist = None
        self.config = config
        self.data = self.read(spots_df)
        self.nS = self.data.shape[0]
        self.unique_gene_names = None
        self._gamma_bar = None
        self._log_gamma_bar = None
        self._gene_id = None
        self._counts_per_gene = None
        [_, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True,
                                                            return_counts=True)

    def __getstate__(self):
        # set here attributes to be excluded from serialisation (pickling)
        # FYI: https://realpython.com/python-pickle-module/
        # Removing _gamma_bar and _log_gamma_bar because they are delayed
        # But even if they werent, I would remove them anyway because they
        # make the pickle file a lot larger!
        attributes = self.__dict__.copy()
        del attributes['_gamma_bar']
        del attributes['_log_gamma_bar']
        return attributes

    # -------- PROPERTIES -------- #
    @property
    def gene_id(self):
        return self._gene_id

    @gene_id.setter
    def gene_id(self, val):
        self._gene_id = val.astype(np.int32)

    @property
    def counts_per_gene(self):
        return self._counts_per_gene

    @counts_per_gene.setter
    def counts_per_gene(self, val):
        self._counts_per_gene = val.astype(np.int32)

    @property
    def gamma_bar(self):
        return self._gamma_bar

    @property
    def log_gamma_bar(self):
        return self._log_gamma_bar

    @property
    def xy_coords(self):
        lst = list(zip(*[self.data.x, self.data.y]))
        return np.array(lst, dtype=np.float32)

    @property
    def parent_cell_prob(self):
        return self._parent_cell_prob

    @parent_cell_prob.setter
    def parent_cell_prob(self, val):
        self._parent_cell_prob = val

    @property
    def parent_cell_id(self):
        return self._parent_cell_id

    @parent_cell_id.setter
    def parent_cell_id(self, val):
        self._parent_cell_id = val

    # -------- METHODS -------- #
    def read(self, spots_df):
        # No need for x_global, y_global to be in the spots_df at first place.
        # Instead of renaming here, you could just use 'x' and 'y' when you
        # created the spots_df
        spots_df = spots_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

        # remove a gene if it is on the exclude list
        exclude_genes = self.config['exclude_genes']
        gene_mask = [True if d not in exclude_genes else False for d in spots_df.target]
        spots_df = spots_df.loc[gene_mask]
        return spots_df.rename_axis('spot_id').rename(columns={'target': 'gene_name'})

    def cells_nearby(self, cells: Cells) -> np.array:
        spotYX = self.data[['y', 'x']].values

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        # maybe return Dist instead of setting the attribute here
        nbrs = cells.nn()
        Dist, neighbors = nbrs.kneighbors(spotYX)
        self.Dist = Dist.astype(np.float32)

        # last column is for misreads.
        neighbors[:, -1] = 0
        return neighbors.astype(np.int32)

    def ini_cellProb(self, neighbors, cfg):
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

    def mvn_loglik(self, data, cell_label, cells):
        centroids = cells.centroid.values[cell_label]
        covs = cells.cov[cell_label]
        # param = list(zip(*[data, centroids, covs]))
        # out = [self.loglik_contr(p) for i, p in enumerate(param)]
        # out_2 = [multivariate_normal.logpdf(p[0], p[1], p[2]) for i, p in enumerate(param)]
        out = self.multiple_logpdfs(data, centroids, covs)
        return out

    def multiple_logpdfs(self, x, means, covs):
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
        Gene counts for the zero expressing class
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

    @dask.delayed
    def gammaExpectation(self, rho, beta):
        """
        :param r:
        :param b:
        :return: Expectetation of a rv X following a Gamma(r,b) distribution with pdf
        f(x;\alpha ,\beta )= \frac{\beta^r}{\Gamma(r)} x^{r-1}e^{-\beta x}
        """
        r = rho[:, :, None]
        return r / beta

    @dask.delayed
    def logGammaExpectation(self, rho, beta):
        r = rho[:, :, None]
        return scipy.special.psi(r) - np.log(beta)


# ----------------------------------------Class: SingleCell--------------------------------------------------- #
class SingleCell(object):
    def __init__(self, scdata: pd.DataFrame, genes: np.array, config):
        self.isMissing = None  # Will be set to False is single cell data are assumed known and given as an input
        # otherwise, if they are unknown, this will be set to True and the algorithm will
        # try to estimate them
        # self.raw_data = self._raw_data(scdata, genes)
        self.config = config
        self._mean_expression, self._log_mean_expression = self._setup(scdata, genes, self.config)

    def _setup(self, scdata, genes, config):
        """
        calcs the mean (and the log-mean) gene counts per cell type. Note that
        some hyperparameter values have been applied before those means are derived.
        These hyperparameters and some bacic cleaning takes part in the functions
        called herein
        """
        if scdata is None:
            datatypes_logger.info('Single Cell data are missing. Cannot determine meam expression per cell class.')
            datatypes_logger.info('We will try to estimate the array instead')
            datatypes_logger.info('Starting point is a diagonal array of size numGenes-by-numGenes')
            # expr = self._naive(scdata, genes)
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
        assert self._mean_expression.columns[-1] == 'Zero', "Last column should be the Zero class"
        return self._mean_expression

    @property
    def log_mean_expression(self):
        assert self._log_mean_expression.columns[-1] == 'Zero', "Last column should be the Zero class"
        return self._log_mean_expression

    @property
    def genes(self):
        return self.mean_expression.index.values

    @property
    def classes(self):
        return self.mean_expression.columns.values

    ## Helper functions ##
    def _set_axes(self, df):
        df = df.rename_axis("class_name", axis="columns").rename_axis('gene_name')
        return df

    def _remove_zero_cols(self, df):
        """
        Removes zero columns (ie if a column is populated by zeros only, then it is removed)
        :param da:
        :return:
        """
        out = df.loc[:, (df != 0).any(axis=0)]
        return out

    def _helper(self, expr):
        # order by column name
        expr = expr.copy().sort_index(axis=0).sort_index(axis=1, key=natsort_keygen(key=lambda y: y.str.lower()))

        # append at the end the Zero class
        expr['Zero'] = np.zeros([expr.shape[0], 1])

        # expr = self.config['Inefficiency'] * arr
        me = expr.rename_axis('gene_name').rename_axis("class_name", axis="columns")  # mean expression
        me = me + self.config['SpotReg']  # add the regularization parameter
        lme = np.log(me)  # log mean expression
        return me, lme

    def _gene_expressions(self, fitted, scale):
        """
        Finds the expected mean gene counts. The prior *IS NOT* taken
        into account. We use data evidence only
        For the zero class only the prior is used *AND NOT* data
        evidence.
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

    def _raw_data(self, scdata, genes):
        """
        Reads the raw single data, filters out any genes outside the gene panel and then it
        groups by the cell type
        """
        assert np.all(scdata >= 0), "Single cell dataframe has negative values"
        datatypes_logger.info(
            'Single cell data passed-in have %d genes and %d cells' % (scdata.shape[0], scdata.shape[1]))

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
        # logger.info('******************************************************')
        # logger.info('*************** DIAGONAL SINGLE CELL DATA ***************')
        # logger.info('******************************************************')
        nG = len(genes)
        mgc = self.config[
            'mean_gene_counts_per_class']  # the avg gene count per cell. Better expose that so it can be set by the user.
        arr = mgc * np.eye(nG)
        labels = ['class_%d' % (i + 1) for i, _ in enumerate(genes)]
        df = pd.DataFrame(arr).set_index(genes)
        df.columns = labels
        return df


# ---------------------------------------- Class: CellType --------------------------------------------------- #
class CellType(object):
    def __init__(self, single_cell, config):
        assert single_cell.classes[-1] == 'Zero', "Last label should be the Zero class"
        self._names = single_cell.classes
        self._alpha = None
        self.config = config
        self.single_cell_data_missing = single_cell.isMissing

    @property
    def names(self):
        assert self._names[-1] == 'Zero', "Last label should be the Zero class"
        return self._names

    @property
    def nK(self):
        return len(self.names)

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        self._alpha = val

    @property
    def pi_bar(self):
        return self.alpha / self.alpha.sum()

    @property
    def logpi_bar(self):
        return scipy.special.psi(self.alpha) - scipy.special.psi(self.alpha.sum())

    @property
    def prior(self):
        return self.pi_bar

    @property
    def log_prior(self):
        if self.single_cell_data_missing or self.config['cell_type_prior'] == 'weighted':
            return self.logpi_bar
        else:
            return np.log(self.prior)

    def size(self, cells):
        """
        calcs the size of a cell class, ie how many members (ie cells) each cell type has
        """
        return cells.classProb.sum(axis=0)

    def ini_prior(self):
        self.alpha = self.ini_alpha()

    def ini_alpha(self):
        ones = np.ones(self.nK - 1)
        return np.append(ones, sum(ones)).astype(np.float32)
