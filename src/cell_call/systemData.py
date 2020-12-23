import numpy as np
import pandas as pd
import xarray as xr
from sklearn.neighbors import NearestNeighbors
import os
import numpy_groupies as npg
import time
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


class Cells(object):
    # Get rid of the properties where not necessary!!
    def __init__(self, config):
        self.config = config
        self.cell_props = read_image_objects(self.config)
        self.classProb = None

    @property
    def yx_coords(self):
        coords = [d for d in zip(self.cell_props.y.values, self.cell_props.x.values) if not np.isnan(d).any()]
        return np.array(coords)

    @property
    def cell_id(self):
        mask = ~np.isnan(self.cell_props.y.values) & ~np.isnan(self.cell_props.x.values)
        return self.cell_props.cell_id.values[mask]

    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.yx_coords.data)
        return nbrs

    def geneCount(self, spots):
        '''
        Produces a matrix numCells-by-numGenes where element at position (c,g) keeps the expected
        number of gene g  in cell c.
        :param spots:
        :return:
        '''
        start = time.time()
        nC = self.yx_coords.shape[0] + 1
        nG = spots.gene_panel.shape[0]
        # cell_id = self.cell_id
        # _id = np.append(cell_id, cell_id.max()+1)
        _id = self.cell_props.index.tolist()
        nN = spots.call.neighbors.shape[1]
        CellGeneCount = np.zeros([nC, nG])

        name = spots.gene_panel.index.values
        ispot = spots.data.gene_id.values
        for n in range(nN - 1):
            c = spots.call.neighbors.loc[:, n].values
            # c = spots.neighboring_cells['id'].sel(neighbor=n).values
            group_idx = np.vstack((c[None, :], ispot[None, :]))
            a = spots.call.cell_prob.loc[:, n]
            accumarray = npg.aggregate(group_idx, a, func="sum", size=(nC, nG))
            CellGeneCount = CellGeneCount + accumarray
        end = time.time()
        print('time in geneCount: ', end - start)
        CellGeneCount = xr.DataArray(CellGeneCount, coords=[_id, name], dims=['cell_id', 'gene_name'])
        # self.CellGeneCount = CellGeneCount
        return CellGeneCount

    def geneCountsPerKlass(self, single_cell_data, egamma, ini):
        temp = self.classProb * self.cell_props.area_factor.to_xarray() * egamma
        temp = temp.sum(dim='cell_id')
        # if you want to calc temp with einsum:
        # temp = np.einsum('ck, c, cgk -> kg', self.classProb, self.cell_props.area_factor.to_xarray(), egamma)
        ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])
        TotPredicted = ClassTotPredicted.drop('Zero', dim='class_name').sum(dim='class_name')
        return TotPredicted

    def geneCountsPerKlass_v2(self, single_cell_data, ini, rho, beta):
        temp = np.einsum('ck, c, cg, cgk -> gk', self.classProb.data, self.cell_props.area_factor.values, rho, 1/beta)
        ClassTotPredicted = np.einsum('gk, gk -> gk', temp, single_cell_data.mean_expression + ini['SpotReg'])
        TotPredicted = np.einsum('gk->g', ClassTotPredicted[:, :-1])  # last class is Zero class. Exclude it from the sum
        # ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])
        # TotPredicted = ClassTotPredicted.drop('Zero', dim='class_name').sum(dim='class_name')
        return TotPredicted


class Cell_prior(object):
    def __init__(self, cell_type):
        # list(dict.fromkeys(cell_type_name))
        self.name = cell_type
        self.nK = self.name.shape[0]
        # Check this....maybe you should divide by K-1
        self.value = np.append([.5 * np.ones(self.nK - 1) / self.nK], 0.5)
        self.logvalue = np.log(self.value)


class Genes(object):
    def __init__(self, spots):
        [gn, ispot, total_spots] = np.unique(spots.data.gene_name.values, return_inverse=True, return_counts=True)
        # self.panel = xr.Dataset({'gene_gamma': xr.DataArray(np.ones(gn.shape), coords=[('gene_name', gn)]),
        #                          'total_spots': xr.DataArray(total_spots, [('gene_name', gn)]),
        #                          },)
        self.panel = pd.DataFrame({'gene_name': gn,
                                   'gene_gamma': 1,  # initial value for gamma is 1
                                   'total_spots': total_spots
                                   })\
            .set_index('gene_name')
        self.ispot = ispot

    def update_gamma(self, cells, spots, single_cell_data, egamma, ini):
        nK = single_cell_data.class_name.shape[0]

        # pSpotZero = spots.zeroKlassProb(klasses, cells)
        TotPredictedZ = spots.TotPredictedZ(self.panel.ispot.data,
                                                cells.classProb.sel({'class_name': 'Zero'}).data)

        TotPredicted = cells.geneCountsPerKlass(single_cell_data, egamma, ini)
        TotPredictedB = np.bincount(spots.geneUniv.ispot.data, spots.neighboring_cells['prob'][:, -1])

        nom = ini['rGene'] + spots.geneUniv.total_spots - TotPredictedB - TotPredictedZ
        denom = ini['rGene'] + TotPredicted
        self.panel.gene_gamma.data = nom / denom


class Spots(object):
    def __init__(self, config):
        self.config = config
        # df = sa.data
        self.data = self.read()
        # df = df.rename_axis('spot_id').rename(columns={'target': 'gene_name'})

        self.call = None
        self._genes = Genes(self)
        self.data['gene_id'] = self._genes.ispot
        self.gene_panel = self._genes.panel
        # self.yxCoords = self.data[['y', 'x']].values

    def read(self):
        saFile = os.path.join(dir_path, self.config['saFile'])

        logger.info('********* Getting spot attributes from %s **********', saFile)
        sa_df = pd.read_csv(saFile)
        if self.config['drop_nan']:
            sa_df = sa_df.dropna()  ##  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE
            logger.info('**** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs****')

        sa_df = sa_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

        # remove a gene if it is on the exclude list
        gene_mask = [True if d not in self.config['exclude_genes'] else False for d in sa_df.target]
        sa_df = sa_df.loc[gene_mask]
        return sa_df.rename_axis('spot_id').rename(columns={'target': 'gene_name'})

    def _neighborCells(self, cells):
        # this needs some clean up.
        spotYX = self.data[['y', 'x']]
        numCells = len(cells.yx_coords)

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotYX)

        # last column is for misreads. Id is dummy id and set to the
        # number of cells (which so-far should always be unallocated)
        neighbors[:, -1] = numCells

        out = pd.DataFrame(neighbors,
                           index=spotYX.index)\
            .rename_axis("neighbour_id", axis="columns")
        return out

    def _cellProb(self, neighbors, cfg):
        nS = self.data.shape[0]
        nN = cfg['nNeighbors'] + 1
        SpotInCell = self.data.label
        assert (np.all(SpotInCell.index == neighbors.index))

        # sanity check (this actually needs to be rewritten)
        mask = np.greater(SpotInCell, 0, where=~np.isnan(SpotInCell))
        sanity_check = neighbors.loc[mask, 0] + 1 == SpotInCell[mask]
        assert ~any(sanity_check), "a spot is in a cell not closest neighbor!"

        pSpotNeighb = np.zeros([nS, nN])
        pSpotNeighb[neighbors.values + 1 == SpotInCell.values[:, None]] = 1
        pSpotNeighb[SpotInCell == 0, -1] = 1

        out = pd.DataFrame(pSpotNeighb,
                           index=neighbors.index,
                           columns=neighbors.columns)

        ## Add a couple of checks here
        return out

    def init_call(self, cells, config):
        neighbors = self._neighborCells(cells)
        my_ds = self._cellProb(neighbors, config)

        out = xr.Dataset({'cell_prob': my_ds, 'neighbors': neighbors})

        assert np.all(out.cell_prob.spot_id.data == out.neighbors.spot_id.data)
        assert np.all(out.cell_prob.neighbour_id.data == out.neighbors.neighbour_id.data)

        self.call = out

    def loglik(self, cells, cfg):
        # meanCellRadius = cells.ds.mean_radius
        mcr = np.mean(np.sqrt(cells.cell_props.area / np.pi)) * 0.5 # This is the meanCellRadius

        # Assume a bivariate normal and calc the likelihood
        # mcr = meanCellRadius.data
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
        spotNeighbours = self.call.neighbors.loc[:, [0, 1, 2]].values

        # get the corresponding probabilities
        neighbourProb = self.call.cell_prob.loc[:, [0, 1, 2]].values

        # prob that a spot belongs to a zero expressing cell
        pSpotZero = np.sum(neighbourProb * pCellZero[spotNeighbours], axis=1)

        # aggregate per gene id
        TotPredictedZ = np.bincount(geneNo, pSpotZero)
        return TotPredictedZ


def read_image_objects(config):
    img_obj = pd.read_csv(config['expanded_cells'])

    meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
    relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

    # append 1 for the misreads
    relCellRadius = np.append(relCellRadius, 1)

    nom = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(config['InsideCellBonus'])) + np.exp(config['InsideCellBonus'])
    denom = np.exp(-0.5) * (1 - np.exp(config['InsideCellBonus'])) + np.exp(config['InsideCellBonus'])
    CellAreaFactor = nom / denom
    areaFactor = CellAreaFactor

    out = pd.DataFrame({'area_factor': areaFactor,
                        'rel_radius': relCellRadius,
                        'area': np.append(img_obj.area, np.nan),
                        # 'mean_radius': meanCellRadius,
                        'x': np.append(img_obj.x.values, np.nan),
                        'y': np.append(img_obj.y.values, np.nan),
                        'cell_id': np.append(img_obj.cell_id.values, img_obj.cell_id.shape[0]+1)})\
        .set_index('cell_id')
    return out


