import os
import json
import numpy as np
import pandas as pd
import numpy_groupies as npg
import scipy
import gc
import numexpr as ne
from sklearn.neighbors import NearestNeighbors
from pciSeq.src.cell_call.utils import read_image_objects
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


class Cells(object):
    # Get rid of the properties where not necessary!!
    def __init__(self, _cells_df, config):
        self.config = config
        self.cell_props = read_image_objects(_cells_df, config)
        self.nC = len(self.cell_props['cell_label'])
        self.classProb = None
        self.class_names = None
        self._prior = None
        self._cov = self.ini_cov()
        self.nu_0 = 15      # need to move that into config.py
        self.rho_1 = self.config['rho_1']    # need to move that into config.py
        self.rho_2 = self.config['rho_2']   # need to move that into config.py
        self._centroid = self.ini_centroids()
        self._gene_counts = None
        self._background_counts = None
        self._alpha = None

    # -------- PROPERTIES -------- #
    @property
    def yx_coords(self):
        coords = [d for d in zip(self.cell_props['y'], self.cell_props['x']) if not np.isnan(d).any()]
        return np.array(coords)

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
        assert val[1:, :].sum() == 0, 'Input array must be zero everywhere apart from the top row'
        self._background_counts = val[0, :]

    @property
    def total_counts(self):
        # tc = self.geneCount.sum(axis=1)
        return self.geneCount.sum(axis=1)

    @property
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, val):
        self._prior = val

    @property
    def log_prior(self):
        return np.log(self.prior)

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        self._alpha = val

    @property
    def centroid(self):
        # lst = list(zip(*[self._centroid['x'], self._centroid['y']]))
        return self._centroid.copy()

    @centroid.setter
    def centroid(self, df):
        self._centroid = df.copy()

    @property
    def cov(self):
        return self._cov

    @cov.setter
    def cov(self, val):
        self._cov = val

    @property
    def sigma_x(self):
        return np.sqrt(self.cov[:, 0, 0])

    @property
    def sigma_y(self):
        return np.sqrt(self.cov[:, 1, 1])

    @property
    def corr(self):
        sigma_x = np.sqrt(self.cov[:, 0, 0])
        sigma_y = np.sqrt(self.cov[:, 1, 1])
        cov_xy = self.cov[:, 0, 1]
        rho = cov_xy / (sigma_x * sigma_y)
        return rho

    @property
    def ellipsoid_attributes(self):
        """
        convenience property that returns a list with
        the ellipoid's centroid, correlation and standard
        deviation across the x-axis and y-axis
        """
        mu = self.centroid.values.tolist()
        # cov = self.cov.tolist()
        sigma_x = self.sigma_x.tolist()
        sigma_y = self.sigma_y.tolist()
        rho = self.corr.tolist()
        # list(zip(self.cells.centroid.values.tolist(), self.cells.sigma_x.tolist(), self.cells.sigma_y.tolist()))
        return list(zip(mu, rho, sigma_x, sigma_y))


    # -------- METHODS -------- #
    def ini_centroids(self):
        d = {'x': self.cell_props['x'], 'y': self.cell_props['y']}
        df = pd.DataFrame(d)
        return df.copy()

    def ini_cov(self):
        mcr = self.dapi_mean_cell_radius()
        cov = mcr * mcr * np.eye(2, 2)
        return np.tile(cov, (self.nC, 1, 1))

    def dapi_mean_cell_radius(self):
        return np.nanmean(np.sqrt(self.cell_props['area'] / np.pi)) * 0.5

    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.yx_coords)
        return nbrs

    def geneCountsPerKlass(self, single_cell_data, egamma, ini):
        temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.alpha, egamma)
        # temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.cell_props['area_factor'], egamma)

        # total counts predicted by all cells of each class (nG, nK)
        ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])

        # total of each gene
        isZero = ClassTotPredicted.columns == 'Zero'
        labels = ClassTotPredicted.columns.values[~isZero]
        TotPredicted = ClassTotPredicted[labels].sum(axis=1)
        return TotPredicted

    def scatter_matrix(self, spots):
        mu_bar = self.centroid.values
        prob = spots.parent_cell_prob[:, :-1]
        id = spots.parent_cell_id[:, :-1]
        xy_spots = spots.xy_coords
        out = self.ini_cov() * self.nu_0
        # out = np.tile(np.eye(2, 2), (self.num_cells, 1, 1))

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
        # the scatter matrix will be a zero squared array. That is fine.
        out[:, 0, 0] = agg_00
        out[:, 1, 1] = agg_11
        out[:, 0, 1] = agg_01
        out[:, 1, 0] = agg_01

        return out

# ----------------------------------------Class: Genes--------------------------------------------------- #
class Genes(object):
    def __init__(self, spots):
        # self.gamma = np.ones(len(spots.unique_gene_names))
        # self.gamma = None
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta = None
        self.nG = len(self.gene_panel)

    @property
    def eta(self):
        return self._eta

    @eta.setter
    def eta(self, val):
        self._eta = val

# ----------------------------------------Class: Spots--------------------------------------------------- #
class Spots(object):
    def __init__(self, spots_df, config):
        self._parent_cell_prob = None
        self._parent_cell_id = None
        self.config = config
        self.data = self.read(spots_df)
        self.nS = self.data.shape[0]
        self.call = None
        self.unique_gene_names = None
        self._gamma_bar = None
        self._log_gamma_bar = None
        [_, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True, return_counts=True)

    # -------- PROPERTIES -------- #
    @property
    def gamma_bar(self):
        return self._gamma_bar.astype(self.config['dtype'])

    @gamma_bar.setter
    def gamma_bar(self, val):
        self._gamma_bar = val.astype(self.config['dtype'])

    @property
    def log_gamma_bar(self):
        return self._log_gamma_bar

    @log_gamma_bar.setter
    def log_gamma_bar(self, val):
        self._log_gamma_bar = val

    @property
    def xy_coords(self):
        lst = list(zip(*[self.data.x, self.data.y]))
        return np.array(lst)

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
        spotYX = self.data[['y', 'x']]

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotYX)

        # last column is for misreads.
        neighbors[:, -1] = 0
        return neighbors

    def ini_cellProb(self, neighbors, cfg):
        nS = self.data.shape[0]
        nN = cfg['nNeighbors'] + 1
        SpotInCell = self.data.label
        # assert (np.all(SpotInCell.index == neighbors.index))

        # sanity check (this actually needs to be rewritten)
        mask = np.greater(SpotInCell, 0, where=~np.isnan(SpotInCell))
        sanity_check = neighbors[mask, 0] + 1 == SpotInCell[mask]
        assert ~any(sanity_check), "a spot is in a cell not closest neighbor!"

        pSpotNeighb = np.zeros([nS, nN])
        pSpotNeighb[neighbors == SpotInCell.values[:, None]] = 1
        pSpotNeighb[SpotInCell == 0, -1] = 1

        ## Add a couple of checks here
        return pSpotNeighb

    def loglik(self, cells, cfg):
        area = cells.cell_props['area'][1:]
        mcr = np.mean(np.sqrt(area / np.pi)) * 0.5  # This is the meanCellRadius

        # Assume a bivariate normal and calc the likelihood
        D = -self.Dist ** 2 / (2 * mcr ** 2) - np.log(2 * np.pi * mcr ** 2)

        # last column (nN-closest) keeps the misreads,
        D[:, -1] = np.log(cfg['MisreadDensity'])

        mask = np.greater(self.data.label, 0, where=~np.isnan(self.data.label))
        D[mask, 0] = D[mask, 0] + cfg['InsideCellBonus']
        return D

    def mvn_loglik(self, data, cell_label, cells):
        centroids = cells.centroid.values[cell_label]
        covs = cells.cov[cell_label]
        param = list(zip(*[data, centroids, covs]))
        out = [self.loglik_contr(p) for i, p in enumerate(param)]
        # out_2 = [multivariate_normal.logpdf(p[0], p[1], p[2]) for i, p in enumerate(param)]
        return out

    def loglik_contr(self, param):
        """
        Contribution of a single datapoint to the bivariate Normal distribution loglikelihood
        Should be the same as multivariate_normal.logpdf(param[0], param[1], param[2])
        """
        x = param[0]
        mu = param[1]
        sigma = param[2]
        k = 2  # dimensionality of a point

        x = x.reshape([2, -1])
        mu = mu.reshape([2, -1])
        assert x.shape == mu.shape == (2, 1), 'Datapoint should have dimension (2, 1)'
        assert sigma.shape == (2, 2), 'Covariance matrix should have shape: (2, 2)'

        sigma_inv = np.linalg.inv(sigma)
        a = - k / 2 * np.log(2 * np.pi)
        b = - 0.5 * np.log(np.linalg.det(sigma))
        c = - 0.5 * np.einsum('ji, jk, ki -> i', x - mu, sigma_inv, x - mu)
        assert len(c) == 1, 'Should be an array with just one element'
        out = a + b + c[0]
        return out

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

    def gammaExpectation(self, rho, beta):
        '''
        :param r:
        :param b:
        :return: Expectetation of a rv X following a Gamma(r,b) distribution with pdf
        f(x;\alpha ,\beta )= \frac{\beta^r}{\Gamma(r)} x^{r-1}e^{-\beta x}
        '''

        # sanity check
        # assert (np.all(rho.coords['cell_id'].data == beta.coords['cell_id'])), 'rho and beta are not aligned'
        # assert (np.all(rho.coords['gene_name'].data == beta.coords['gene_name'])), 'rho and beta are not aligned'

        dtype = self.config['dtype']
        r = rho[:, :, None]
        if dtype == np.float64:
            gamma = np.empty(beta.shape)
            ne.evaluate('r/beta', out=gamma)
            return gamma
        else:
            return (r/beta).astype(dtype)

    def logGammaExpectation(self, rho, beta):
        dtype = self.config['dtype']
        r = rho[:, :, None].astype(dtype)
        if dtype == np.float64:
            logb = np.empty(beta.shape)
            ne.evaluate("log(beta)", out=logb)
            return scipy.special.psi(r) - logb
        else:
            return scipy.special.psi(r) - np.log(beta).astype(dtype)


# ----------------------------------------Class: SingleCell--------------------------------------------------- #
class SingleCell(object):
    def __init__(self, scdata: pd.DataFrame, genes: np.array, config):
        self.config = config
        self._mean_expression, self._log_mean_expression = self._setup(scdata, genes, self.config)

    def _setup(self, scdata, genes, config):
        assert np.all(scdata >= 0), "Single cell dataframe has negative values"
        logger.info(' Single cell data passed-in have %d genes and %d cells' % (scdata.shape[0], scdata.shape[1]))

        logger.info(' Single cell data: Keeping counts for the gene panel of %d only' % len(genes))
        df = scdata.loc[genes]

        # set the axes labels
        df = self._set_axes(df)

        df = self._remove_zero_cols(df.copy())
        dfT = df.T

        logger.info(' Single cell data: Grouping gene counts by cell type. Aggregating function is the mean.')
        expr = dfT.groupby(dfT.index.values).agg('mean').T
        expr['Zero'] = np.zeros([expr.shape[0], 1])
        expr = expr.sort_index(axis=0).sort_index(axis=1)
        expr = config['Inefficiency'] * expr
        me = expr.rename_axis('gene_name').rename_axis("class_name", axis="columns")  # mean expression
        lme = np.log(me + config['SpotReg'])  # log mean expression

        logger.info(' Grouped single cell data have %d genes and %d cell types' % (me.shape[0], me.shape[1]))
        dtype = self.config['dtype']
        return me.astype(dtype), lme.astype(dtype)

    # -------- PROPERTIES -------- #
    @property
    def mean_expression(self):
        return self._mean_expression

    @property
    def log_mean_expression(self):
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





