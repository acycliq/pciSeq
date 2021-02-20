import os
import json
import numpy as np
import pandas as pd
import numpy_groupies as npg
from sklearn.neighbors import NearestNeighbors
from pciSeq.src.cell_call.utils import read_image_objects
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
        self.num_cells = len(self.cell_props['cell_label'])
        self.classProb = None
        self.class_names = None
        self.log_prior = None
        # initialise covariance matrix for all cells
        self._cov = np.tile(8.867 * 8.867 * np.eye(2, 2), (self.num_cells, 1, 1))
        # initialise centroids from the cell segmentation
        self._centroid = {'x': self.cell_props['x'], 'y': self.cell_props['y']}
        self._gene_counts = None

    @property
    def yx_coords(self):
        coords = [d for d in zip(self.cell_props['y'], self.cell_props['x']) if not np.isnan(d).any()]
        return np.array(coords)

    @property
    def centroid(self):
        lst = list(zip(*[self._centroid['x'], self._centroid['y']]))
        return np.array(lst)

    @centroid.setter
    def centroid(self, val):
        # Add a check to prevent cell at position 0 from being mutated
        try:
            centre, cell_id = val
        except ValueError:
            mask = np.arange(self.num_cells)
        else:
            mask = cell_id
        self._centroid['x'][mask] = centre[0]
        self._centroid['y'][mask] = centre[1]

    @property
    def cov(self):
        return self._cov

    @cov.setter
    def cov(self, val):
        try:
            cov, cell_id = val
        except ValueError:
            mask = np.arange(self.num_cells)
        else:
            mask = cell_id
        if cell_id:
            self._cov[mask] = cov



    # @property
    # def cell_id(self):
    #     mask = ~np.isnan(self.cell_props.y) & ~np.isnan(self.cell_props.x)
    #     return self.cell_props.cell_id[mask]

    def prior(self, cell_type):
        name = cell_type
        nK = name.shape[0]
        # Check this....maybe you should divide my K-1
        value = np.append([.5 * np.ones(nK - 1) / nK], 0.5)
        return np.log(value)

    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.yx_coords)
        return nbrs

    def geneCount(self, spots):
        '''
        Produces a matrix numCells-by-numGenes where element at position (c,g) keeps the expected
        number of gene g  in cell c.
        :param spots:
        :return:
        '''
        start = time.time()
        nC = self.num_cells
        nG = len(spots.unique_gene_names)
        # cell_id = self.cell_id
        # _id = np.append(cell_id, cell_id.max()+1)
        # _id = self.cell_props['cell_id']
        nN = self.config['nNeighbors'] + 1
        CellGeneCount = np.zeros([nC, nG])

        # name = spots.gene_panel.index.values
        spot_id = spots.gene_id
        for n in range(nN - 1):
            c = spots.adj_cell_id[:, n]
            # c = spots.neighboring_cells['id'].sel(neighbor=n).values
            group_idx = np.vstack((c[None, :], spot_id[None, :]))
            a = spots.adj_cell_prob[:, n]
            accumarray = npg.aggregate(group_idx, a, func="sum", size=(nC, nG))
            CellGeneCount = CellGeneCount + accumarray
        end = time.time()
        # print('time in geneCount: ', end - start)
        # CellGeneCount = xr.DataArray(CellGeneCount, coords=[_id, name], dims=['cell_id', 'gene_name'])
        # self.CellGeneCount = CellGeneCount
        return CellGeneCount

    def geneCountsPerKlass(self, single_cell_data, egamma, ini):
        # temp = self.classProb * self.cell_props.area_factor.to_xarray() * egamma
        # temp = temp.sum(dim='cell_id')
        # if you want to calc temp with einsum:
        temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.cell_props['area_factor'], egamma)
        ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])
        TotPredicted = ClassTotPredicted.drop('Zero', dim='class_name').sum(dim='class_name')
        return TotPredicted


class Genes(object):
    def __init__(self, spots):
        self.gamma = np.ones(len(spots.unique_gene_names))
        self.gene_names = spots.unique_gene_names

    def update_gamma(self, cells, spots, single_cell_data, egamma, ini):
        TotPredictedZ = spots.TotPredictedZ(self.panel.spot_id.data,
                                                cells.classProb.sel({'class_name': 'Zero'}).data)

        TotPredicted = cells.geneCountsPerKlass(single_cell_data, egamma, ini)
        TotPredictedB = np.bincount(spots.geneUniv.spot_id.data, spots.neighboring_cells['prob'][:, -1])

        nom = ini['rGene'] + spots.geneUniv.total_spots - TotPredictedB - TotPredictedZ
        denom = ini['rGene'] + TotPredicted
        self.panel.gene_gamma.data = nom / denom


class Spots(object):
    def __init__(self, spots_df, config):
        self.config = config
        self.data = self.read(spots_df)
        self.call = None
        self.adj_cell_prob = None
        self.adj_cell_id = None
        self.unique_gene_names = None
        self.gene_id = None
        self.counts_per_gene = None
        self._unique_genes()
        # self._genes = Genes(self)
        # self.data['gene_id'] = self._genes.spot_id
        # self.gene_panel = self._genes.panel

    @property
    def xy_coords(self):
        lst = list(zip(*[self.data.x, self.data.y]))
        return np.array(lst)

    def _unique_genes(self):
        [self.unique_gene_names, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True, return_counts=True)
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

    def _neighborCells(self, cells):
        # this needs some clean up.
        spotYX = self.data[['y', 'x']]
        # numCells = cells.num_cells

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotYX)

        # last column is for misreads.
        neighbors[:, -1] = 0

        return neighbors

    def _cellProb(self, neighbors, cfg):
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

    def init_call(self, cells, config):
        self.adj_cell_id = self._neighborCells(cells)
        self.adj_cell_prob = self._cellProb(self.adj_cell_id, config)

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
        centroids = cells.centroid[cell_label]
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
        assert x.shape == mu.shape == (2,1), 'Datapoint should have dimension (2, 1)'
        assert sigma.shape == (2, 2), 'Covariance matrix should have shape: (2, 2)'

        sigma_inv = np.linalg.inv(sigma)
        a = - k/2 * np.log(2 * np.pi)
        b = - 0.5 * np.log(np.linalg.det(sigma))
        c = - 0.5 * np.einsum('ji, jk, ki -> i', x-mu, sigma_inv, x-mu)
        assert len(c) == 1, 'Should be an array with just one element'
        out = a + b + c[0]
        return out

    def TotPredictedZ(self, geneNo, pCellZero):
        '''
        ' given a vector
        :param spots:
        :return:
        '''

        # for each spot get the ids of the 3 nearest cells
        spotNeighbours = self.adj_cell_id[:, :-1]

        # get the corresponding probabilities
        neighbourProb = self.adj_cell_prob[:, :-1]

        # prob that a spot belongs to a zero expressing cell
        pSpotZero = np.sum(neighbourProb * pCellZero[spotNeighbours], axis=1)

        # aggregate per gene id
        TotPredictedZ = np.bincount(geneNo, pSpotZero)
        return TotPredictedZ




