import numpy as np
import pandas as pd
import xarray as xr
from sklearn.neighbors import NearestNeighbors
import os
import numpy_groupies as npg
from src.cell_call.utils import read_image_objects
import time
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


class Cells(object):
    # Get rid of the properties where not necessary!!
    def __init__(self, config):
        self.config = config
        self.cell_props = read_image_objects(self.config)
        self.num_cells = len(self.cell_props['cell_id'])
        self.classProb = None
        self.className = None

    @property
    def yx_coords(self):
        coords = [d for d in zip(self.cell_props['y'], self.cell_props['x']) if not np.isnan(d).any()]
        return np.array(coords)

    # @property
    # def cell_id(self):
    #     mask = ~np.isnan(self.cell_props.y) & ~np.isnan(self.cell_props.x)
    #     return self.cell_props.cell_id[mask]

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
        print('time in geneCount: ', end - start)
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


class Cell_prior(object):
    def __init__(self, cell_type):
        # list(dict.fromkeys(cell_type_name))
        self.name = cell_type
        self.nK = self.name.shape[0]
        # Check this....maybe you should divide my K-1
        self.value = np.append([.5 * np.ones(self.nK - 1) / self.nK], 0.5)
        self.logvalue = np.log(self.value)


class Genes(object):
    def __init__(self, spots):
        # [gn, spot_id, total_spots] = np.unique(spots.data.gene_name.values, return_inverse=True, return_counts=True)
        self.gamma = np.ones(len(spots.unique_gene_names))
        self.gene_names = spots.unique_gene_names


        # self.panel = pd.DataFrame({'gene_name': gn,
        #                            'gene_gamma': 1,  # initial value for gamma is 1
        #                            'total_spots': total_spots
        #                            })\
        #     .set_index('gene_name')
        # # self.spot_id = spot_id

    def update_gamma(self, cells, spots, single_cell_data, egamma, ini):
        TotPredictedZ = spots.TotPredictedZ(self.panel.spot_id.data,
                                                cells.classProb.sel({'class_name': 'Zero'}).data)

        TotPredicted = cells.geneCountsPerKlass(single_cell_data, egamma, ini)
        TotPredictedB = np.bincount(spots.geneUniv.spot_id.data, spots.neighboring_cells['prob'][:, -1])

        nom = ini['rGene'] + spots.geneUniv.total_spots - TotPredictedB - TotPredictedZ
        denom = ini['rGene'] + TotPredicted
        self.panel.gene_gamma.data = nom / denom


class Spots(object):
    def __init__(self, config):
        self.config = config
        self.data = self.read()
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

    def _unique_genes(self):
        [self.unique_gene_names, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True, return_counts=True)
        return self.data.gene_name.values

    def read(self):
        spotsFile = os.path.join(dir_path, self.config['spotsFile'])

        logger.info('********* Getting spot attributes from %s **********', spotsFile)
        spots_df = pd.read_csv(spotsFile)
        if 'drop_nan' in self.config.keys() and self.config['drop_nan']:
            spots_df = spots_df.dropna()  ##  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE
            logger.info('**** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs****')

        spots_df = spots_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

        # remove a gene if it is on the exclude list
        gene_mask = [True if d not in self.config['exclude_genes'] else False for d in spots_df.target]
        spots_df = spots_df.loc[gene_mask]
        return spots_df.rename_axis('spot_id').rename(columns={'target': 'gene_name'})

    def _neighborCells(self, cells):
        # this needs some clean up.
        spotYX = self.data[['y', 'x']]
        numCells = cells.num_cells

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotYX)

        # last column is for misreads. Id is dummy id and set to the
        # number of cells (which so-far should always be unallocated)
        neighbors[:, -1] = numCells

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
        pSpotNeighb[neighbors + 1 == SpotInCell.values[:, None]] = 1  # neighbors is 0-based whereas SpotInCell 1-based
        pSpotNeighb[SpotInCell == 0, -1] = 1

        ## Add a couple of checks here
        return pSpotNeighb

    def init_call(self, cells, config):
        self.adj_cell_id = self._neighborCells(cells)
        self.adj_cell_prob = self._cellProb(self.adj_cell_id, config)

    def loglik(self, cells, cfg):
        # meanCellRadius = cells.ds.mean_radius
        area = cells.cell_props['area'][:-1]
        mcr = np.mean(np.sqrt(area / np.pi)) * 0.5  # This is the meanCellRadius

        # Assume a bivariate normal and calc the likelihood
        D = -self.Dist ** 2 / (2 * mcr ** 2) - np.log(2 * np.pi * mcr ** 2)

        # last column (nN-closest) keeps the misreads,
        D[:, -1] = np.log(cfg['MisreadDensity'])

        mask = np.greater(self.data.label, 0, where=~np.isnan(self.data.label))
        D[mask, 0] = D[mask, 0] + cfg['InsideCellBonus']
        print('in loglik')
        return D

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




