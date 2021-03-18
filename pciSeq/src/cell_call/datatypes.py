import os
import json
import numpy as np
import pandas as pd
import numpy_groupies as npg
from sklearn.neighbors import NearestNeighbors
from pciSeq.src.cell_call.utils import read_image_objects
from typing import Tuple
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
from scipy.stats import multivariate_normal
import time
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
        # initialise covariance matrix for all cells
        self._cov = self.ini_cov()
        self.nu_0 = 15      # need to move that into config.py
        self.rho_1 = 100    # need to move that into config.py
        self.rho_2 = 100   # need to move that into config.py
        # self._cov = np.tile(8.867 * 8.867 * np.eye(2, 2), (self.num_cells, 1, 1))
        # initialise centroids from the cell segmentation
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
        # return np.array(list(zip(rho, sigma_x, sigma_y)))

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
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        self._alpha = val

    @property
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, val):
        self._prior = val

    @property
    def log_prior(self):
        return np.log(self.prior)

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

    def refresh_me(self, spots) -> None:
        self.geneCount_upd(spots)
        self.centroid_upd(spots)
        self.cov_upd(spots)
        # self.corr
        # logger.info('refreshed!')

    # def centroid_upd(self, spots):
    #     # get the total gene counts per cell
    #     N_c = self.total_counts
    #
    #     xy_spots = spots.xy_coords
    #     prob = spots.parent_cell_prob
    #     n = self.config['nNeighbors'] + 1
    #
    #     # mulitply the x coord of the spots by the cell prob
    #     a = np.tile(xy_spots[:, 0], (n, 1)).T * prob
    #
    #     # mulitply the y coord of the spots by the cell prob
    #     b = np.tile(xy_spots[:, 1], (n, 1)).T * prob
    #
    #     # aggregated x and y coordinate
    #     idx = spots.parent_cell_id
    #     _x_agg = npg.aggregate(idx.ravel(), a.ravel())
    #     _y_agg = npg.aggregate(idx.ravel(), b.ravel())
    #
    #     x_agg = np.zeros(N_c.shape)
    #     mask = np.arange(len(_x_agg))
    #     x_agg[mask] = _x_agg
    #
    #     y_agg = np.zeros(N_c.shape)
    #     mask = np.arange(len(_y_agg))
    #     y_agg[mask] = _y_agg
    #
    #     # get the estimated cell centers
    #     x_bar = np.nan * np.ones(N_c.shape)
    #     y_bar = np.nan * np.ones(N_c.shape)
    #     x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
    #     y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]
    #     # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
    #     xy_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T)))
    #
    #     # if you have a value for the estimated centroid use that, otherwise
    #     # use the initial (starting values) centroids
    #     ini_cent = self.ini_centroids()
    #     xy_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y']])))
    #
    #     # # sanity check. NaNs or Infs should appear together
    #     # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
    #     # use the fitted centroids where possible otherwise use the initial ones
    #     xy_bar[np.isfinite(x_bar)] = xy_bar_fitted[np.isfinite(x_bar)]
    #     self.centroid = pd.DataFrame(xy_bar, columns=['x', 'y'])
    #     # print(np.array(list(zip(x_bar.T, y_bar.T))))
    #
    # def cov_upd(self, spots):
    #     # first get the scatter matrix
    #     S = self.scatter_matrix(spots)  # sample sum of squares
    #     cov_0 = self.ini_cov()
    #     nu_0 = self.nu_0
    #     S_0 = cov_0 * nu_0 # prior sum of squarea
    #     N_c = self.total_counts
    #     d = 2
    #     denom = np.ones([self.num_cells, 1, 1])
    #     denom[:, 0, 0] = N_c + nu_0
    #     # Note: need to add code to handle the case N_c + nu_0 <= d + 2
    #     cov = (S + S_0) / denom
    #     self.cov = cov

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

        # # The following is for the Ledoit Wolf shrinkage. It implements pi from
        # # page 13 of http://www.ledoit.net/honey.pdf
        # denom = np.ones([self.nC, 1, 1])
        # denom[:, 0, 0] = self.total_counts
        # sample_cov = np.tile(np.zeros([2, 2]), (self.nC, 1, 1))
        # mask = self.total_counts > 0
        # sample_cov[mask] = out[mask]/denom[mask]
        #
        #
        # sq_el_00 = prob * x_centered**2 * x_centered**2  # contribution to the scatter matrix's [0, 0] element
        # sq_el_11 = prob * y_centered**2 * y_centered**2  # contribution to the scatter matrix's off-diagonal element
        # sq_el_01 = prob * x_centered**2 * y_centered**2  # contribution to the scatter matrix's [1, 1] element
        #
        # # Aggregate all contributions to get the scatter matrix
        # sq_agg_00 = npg.aggregate(id.ravel(), sq_el_00.ravel(), size=self.nC)
        # sq_agg_11 = npg.aggregate(id.ravel(), sq_el_11.ravel(), size=self.nC)
        # sq_agg_01 = npg.aggregate(id.ravel(), sq_el_01.ravel(), size=self.nC)
        #
        # sq_sample_cov = np.zeros(out.shape)
        # sq_sample_cov[mask, 0, 0] = sq_agg_00[mask]/self.total_counts[mask]
        # sq_sample_cov[mask, 1, 1] = sq_agg_11[mask]/self.total_counts[mask]
        # sq_sample_cov[mask, 0, 1] = sq_agg_01[mask]/self.total_counts[mask]
        # sq_sample_cov[mask, 1, 0] = sq_agg_01[mask]/self.total_counts[mask]
        #
        # pi_mat = sq_sample_cov - sample_cov**2
        # #pi hat
        # pi_hat = [d.sum() for d in pi_mat]
        #
        # #rho hat
        # rho_hat = [d[0, 0] + d[1, 1] for d in pi_mat]
        #
        # # gamma-hat
        # prior = self.ini_cov()
        # gamma_hat = [np.linalg.norm(sample_cov[i] - prior[i], "fro") ** 2 for i in range(self.nC)]
        #
        # # shrinkage constant
        # kappa = (np.array(pi_hat) - np.array(rho_hat)) / np.array(gamma_hat)
        #
        # shrink = np.zeros(self.nC)
        # shrink[mask] = np.maximum(0, np.minimum(1, kappa[mask] / self.total_counts[mask]))
        return out

    def ledoit_wolf(self, spots, sample_cov):
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

        sq_el_00 = prob * x_centered**2 * x_centered**2  # contribution to the scatter matrix's [0, 0] element
        sq_el_11 = prob * y_centered**2 * y_centered**2  # contribution to the scatter matrix's off-diagonal element
        sq_el_01 = prob * x_centered**2 * y_centered**2  # contribution to the scatter matrix's [1, 1] element

        # Aggregate all contributions to get the scatter matrix
        sq_agg_00 = npg.aggregate(id.ravel(), sq_el_00.ravel(), size=self.nC)
        sq_agg_11 = npg.aggregate(id.ravel(), sq_el_11.ravel(), size=self.nC)
        sq_agg_01 = npg.aggregate(id.ravel(), sq_el_01.ravel(), size=self.nC)

        cov_0 = self.ini_cov()
        nu_0 = self.nu_0
        S_0 = cov_0 * nu_0  # prior sum of squarea

        sq_scatter_mat = np.zeros(out.shape)
        sq_scatter_mat[:, 0, 0] = sq_agg_00
        sq_scatter_mat[:, 1, 1] = sq_agg_11
        sq_scatter_mat[:, 0, 1] = sq_agg_01
        sq_scatter_mat[:, 1, 0] = sq_agg_01

        sq_scatter_mat = sq_scatter_mat + S_0
        denom = np.ones([self.nC, 1, 1])
        denom[:, 0, 0] = self.total_counts + nu_0
        sq_sample_cov = sq_scatter_mat / denom

        pi_mat = sq_sample_cov - sample_cov ** 2
        # pi hat
        pi_hat = np.array([d.sum() for d in pi_mat])

        # rho hat
        rho_hat = np.array([d[0, 0] + d[1, 1] for d in pi_mat])

        # gamma-hat
        prior = self.ini_cov()
        gamma_hat = np.array([np.linalg.norm(sample_cov[i] - prior[i], "fro") ** 2 for i in range(self.nC)])

        # shrinkage constant
        mask = self.total_counts > 0
        kappa = np.zeros(self.nC)
        kappa[mask] = (pi_hat[mask] - rho_hat[mask]) / gamma_hat[mask]

        shrink = np.ones(self.nC) * 0.5
        shrink[mask] = np.maximum(0, np.minimum(1, kappa[mask] / self.total_counts[mask]))
        return shrink

    def stein(self, x: np.array) -> np.array:
        """
        x: The covariance matrix
        """
        # stein_1 = self.stein_test()
        _shape = x.shape
        assert _shape[1] == _shape[2] == 2, "Input must be array-like of 2-by-2 arrays"
        p = _shape[1]

        [eigval, eigvec] = np.linalg.eig(x)
        # D = np.eye(p)
        # np.fill_diagonal(D, eigval)
        # sigma = eigvec @ D @ eigvec.transpose()

        # p = 2
        # n = 15
        # stein = np.zeros([p, p])

        # Equation 8 from https://link.springer.com/content/pdf/10.1007/s00180-016-0672-4.pdf
        N_c = self.total_counts
        mask = eigval[:, 0] != eigval[:, 1]
        a_0 = np.zeros(self.nC)
        a_1 = np.zeros(self.nC)
        a_0[mask] = N_c[mask] - p + 1 + 2 * eigval[mask, 0] * 1 / (eigval[mask, 0] - eigval[mask, 1])
        a_1[mask] = N_c[mask] - p + 1 + 2 * eigval[mask, 1] * 1 / (eigval[mask, 1] - eigval[mask, 0])

        # Equation 9 from https://link.springer.com/content/pdf/10.1007/s00180-016-0672-4.pdf
        phi = np.zeros(x.shape)
        phi[:, 0, 0] = eigval[:, 0]
        phi[:, 1, 1] = eigval[:, 1]
        phi[mask, 0, 0] = eigval[mask, 0] / a_0[mask]
        phi[mask, 1, 1] = eigval[mask, 1] / a_1[mask]

        # the following is simply sigma = eigvec @ phi @ eigvec.transpose()
        sigma = np.einsum('cpn, cnm, cqm -> cpq', eigvec, phi, eigvec)
        return sigma

    def stein_test(self):
        cov = [[65.5299, -22.4162], [-22.4162, 80.3575]]
        p = 2
        N_c = 14
        [eigval, eigvec] = np.linalg.eig(cov)
        D = np.eye(p)
        np.fill_diagonal(D, eigval)
        sigma = eigvec @ D @ eigvec.transpose()
        phi = np.zeros([2,2])
        a_0 = N_c - p + 1 + 2 * eigval[0] * 1 / (eigval[0] - eigval[1])
        a_1 = N_c - p + 1 + 2 * eigval[1] * 1 / (eigval[1] - eigval[0])

        phi[0, 0] = eigval[0] / a_0
        phi[1, 1] = eigval[1] / a_1

        stein = eigvec @ phi @ eigvec.transpose()
        return stein




    # def stein(self. spots):


    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.yx_coords)
        return nbrs

    def geneCountsPerKlass(self, single_cell_data, egamma, ini):
        # temp = self.classProb * self.cell_props.area_factor.to_xarray() * egamma
        # temp = temp.sum(dim='cell_id')
        # if you want to calc temp with einsum:
        temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.alpha, egamma)
        # temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.cell_props['area_factor'], egamma)

        # total counts predicted by all cells of each class (nG, nK)
        ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])

        # total of each gene
        TotPredicted = ClassTotPredicted.drop('Zero', dim='class_name').sum(dim='class_name')
        return TotPredicted


class Genes(object):
    def __init__(self, spots):
        # self.gamma = np.ones(len(spots.unique_gene_names))
        # self.gamma = None
        self._eta = None
        self.gene_names = spots.unique_gene_names
        self.nG = len(self.gene_names)

    @property
    def eta(self):
        return self._eta

    @eta.setter
    def eta(self, val):
        self._eta = val


class Spots(object):
    def __init__(self, spots_df, config):
        self._parent_cell_prob = None
        self._parent_cell_id = None
        self.config = config
        self.data = self.read(spots_df)
        self.nS = self.data.shape[0]
        self.call = None
        self.unique_gene_names = None
        self.gene_id = None
        self.counts_per_gene = None
        self._unique_genes()
        self._gamma_bar = None
        self._log_gamma_bar = None

    @property
    def gamma_bar(self):
        return self._gamma_bar

    @gamma_bar.setter
    def gamma_bar(self, val):
        self._gamma_bar = val

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

    # def update_cell_prob(self, new_assignments, cell_obj):
    #     # Updates the parent cell probabilities
    #     self._adj_cell_prob = new_assignments
    #     cell_obj.refresh_me(self)  # Since the spot-to-cell prob changed you have to change the cell gene counts too

    def _unique_genes(self):
        [self.unique_gene_names, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values,
                                                                                 return_inverse=True,
                                                                                 return_counts=True)
        return self.data.gene_name.values

    def read(self, spots_df):
        # tempdir = self.config['PREPROCESS']['temp']
        # spotsFile = os.path.join(tempdir, '_spots.csv')

        # logger.info('********* Getting spot attributes from %s **********', spotsFile)
        # spots_df = pd.read_csv(spotsFile)
        # spots_df = spots_df.sample(frac=1).reset_index(drop=True)

        if 'drop_nan' in self.config.keys() and self.config['drop_nan']:
            spots_df = spots_df.dropna()  ##  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE
            logger.info('**** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs****')

        spots_df = spots_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

        # remove a gene if it is on the exclude list
        exclude_genes = self.config['exclude_genes']
        gene_mask = [True if d not in exclude_genes else False for d in spots_df.target]
        spots_df = spots_df.loc[gene_mask]
        return spots_df.rename_axis('spot_id').rename(columns={'target': 'gene_name'})

    def cells_nearby(self, cells: Cells) -> np.array:
        # this needs some clean up.
        spotYX = self.data[['y', 'x']]
        # numCells = cells.num_cells

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotYX)

        # last column is for misreads.
        neighbors[:, -1] = 0

        return neighbors

    def ini_cellProb(self, neighbors, cfg):
        # Note: Something is not right here.
        # The total sum of the return value should be the same as the number od spots but it is NOT!!
        # It doesnt seem to be a very crucial bug though.
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

    # def init_call(self, cells, config):
    #     self.adj_cell_id = self.cells_nearby(cells)
    #     ini_prob = self.ini_cellProb(self.adj_cell_id, config)
    #     self.update_cell_prob(ini_prob, cells)
    #     logger.info('ok')
    #     # self.adj_cell_prob = self.ini_cellProb(self.adj_cell_id, config)

    def loglik(self, cells, cfg):
        # meanCellRadius = cells.ds.mean_radius
        area = cells.cell_props['area'][1:]
        mcr = np.mean(np.sqrt(area / np.pi)) * 0.5  # This is the meanCellRadius

        # Assume a bivariate normal and calc the likelihood
        D = -self.Dist ** 2 / (2 * mcr ** 2) - np.log(2 * np.pi * mcr ** 2)

        # last column (nN-closest) keeps the misreads,
        D[:, -1] = np.log(cfg['MisreadDensity'])

        mask = np.greater(self.data.label, 0, where=~np.isnan(self.data.label))
        D[mask, 0] = D[mask, 0] + cfg['InsideCellBonus']
        # print('in loglik')
        return D

    # def mvn_loglik_par(self, data, mu, sigma):
    #     n = max(1, cpu_count() - 1)
    #     pool = ThreadPool(n)
    #     param = list(zip(*[data, mu, sigma]))
    #     results = pool.map(self.loglik_contr, param)
    #     # out = [self.loglik_contr(data[i], mu[i], sigma[i]) for i, d in enumerate(data)]
    #     pool.close()
    #     pool.join()
    #     return results

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
        '''
        ' given a vector
        :param spots:
        :return:
        '''

        # for each spot get the ids of the 3 nearest cells
        spotNeighbours = self.parent_cell_id[:, :-1]

        # get the corresponding probabilities
        neighbourProb = self.parent_cell_prob[:, :-1]

        # prob that a spot belongs to a zero expressing cell
        pSpotZero = np.sum(neighbourProb * pCellZero[spotNeighbours], axis=1)

        # aggregate per gene id
        TotPredictedZ = np.bincount(geneNo, pSpotZero)
        return TotPredictedZ
