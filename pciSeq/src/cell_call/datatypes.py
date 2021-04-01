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
        self._cov = self.ini_cov()
        self._gene_counts = None
        self._background_counts = None

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

    # -------- METHODS -------- #
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
        # temp = self.classProb * self.cell_props.area_factor.to_xarray() * egamma
        # temp = temp.sum(dim='cell_id')
        # if you want to calc temp with einsum:
        # temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.alpha, egamma)
        temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.cell_props['area_factor'], egamma)

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
