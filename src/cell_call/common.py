import os
import src.cell_call.utils as utils
import numpy as np
import pandas as pd
import xarray as xr
import gc
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


def expected_gamma(cells, spots, ds, ini):
    scaled_mean = cells.ds.area_factor.to_xarray() * ds.mean_expression
    rho = ini['rSpot'] + cells.geneCount(spots)
    beta = ini['rSpot'] + scaled_mean

    expected_gamma = utils.gammaExpectation(rho, beta)
    expected_loggamma = utils.logGammaExpectation(rho, beta)

    del rho
    del beta
    gc.collect()
    del gc.garbage[:]

    return expected_gamma, expected_loggamma


def celltype_assignment(cells, spots, prior, ds, cfg):
    '''
    return a an array of size numCells-by-numCellTypes where element in position [i,j]
    keeps the probability that cell i has cell type j
    :param spots:
    :param config:
    :return:
    '''

    gene_gamma = spots.gene_panel.gene_gamma
    ScaledExp = cells.ds.area_factor.to_xarray() * gene_gamma.to_xarray() * ds.mean_expression + cfg['SpotReg']
    pNegBin = ScaledExp / (cfg['rSpot'] + ScaledExp)
    # contr = utils.nb_negBinLoglik(CellGeneCount[:,:,None], ini['rSpot'], pNegBin)
    cgc = cells.geneCount(spots)
    # x * log(p) + r * log(1 - p)
    contr = utils.negBinLoglik(cgc, cfg['rSpot'], pNegBin)
    # assert np.all(nb_contr == contr)
    wCellClass = np.sum(contr, axis=1) + prior.logvalue
    # pCellClass_old = utils.softmax(wCellClass)
    pCellClass = xr.apply_ufunc(utils.softmax_nolan, wCellClass, 1, 1)

    cells.classProb = pCellClass
    logger.info('Cell 0 is classified as %s with prob %4.8f' % (
        prior.name[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
    logger.info('cell ---> klass probabilities updated')
    return pCellClass


def call_spots(spots, cells, single_cell_data, prior, elgamma, cfg):
    # spot to cell assignment
    nN = spots.call.neighbors.shape[1]
    nS = spots.data.gene_name.shape[0]
    nK = prior.nK
    aSpotCell = np.zeros([nS, nN])
    gn = spots.data.gene_name.values
    expected_spot_counts = single_cell_data['log_mean_expression'].sel({'gene_name': gn}).data
    for n in range(nN - 1):
        spots_name = spots.data.gene_name.values
        single_cell_data['log_mean_expression'].sel({'gene_name': spots_name})

        # get the spots' nth-closest cell
        sn = spots.call.neighbors.loc[:, n].values

        # get the respective cell type probabilities
        cp = cells.classProb.sel({'cell_id': sn}).data

        # multiply and sum over cells
        term_1 = (expected_spot_counts * cp).sum(axis=1)

        # logger.info('genes.spotNo should be something line spots.geneNo instead!!')
        # expectedLog = elgamma.data[spots.call.neighbors[:, n], spots.data.gene_id]
        expectedLog = utils.bi2(elgamma.data, [nS, nK], sn[:, None], spots.data.gene_id.values[:,None])
        term_2 = np.sum(cp * expectedLog, axis=1)
        aSpotCell[:, n] = term_1 + term_2
    wSpotCell = aSpotCell + spots.loglik(cells, cfg)

    # update the prob a spot belongs to a neighboring cell
    pSpotNeighb = utils.softmax2(wSpotCell)
    spots.call.cell_prob.data = pSpotNeighb
    logger.info('spot ---> cell probabilities updated')


def updateGamma(cells, spots, single_cell_data, egamma, ini):
    # Maybe I should rename that to eta (not gamma). In the paper the symbol is eta
    nK = single_cell_data.class_name.shape[0]

    # pSpotZero = spots.zeroKlassProb(klasses, cells)
    TotPredictedZ = spots.TotPredictedZ(spots.data.gene_id.values, cells.classProb.sel({'class_name': 'Zero'}).data)

    TotPredicted = cells.geneCountsPerKlass(single_cell_data, egamma, ini)

    TotPredictedB = np.bincount(spots.data.gene_id.values, spots.call.cell_prob.loc[:, 3].data)

    nom = ini['rGene'] + spots.gene_panel.total_spots - TotPredictedB - TotPredictedZ
    denom = ini['rGene'] + TotPredicted
    res = nom / denom

    # Finally, update gene_gamma
    spots.gene_panel.gene_gamma[res.index] = res.values
    # cells.expectedGamma = nom / denom


def geneCountsPerKlass(cells, single_cell_data, egamma, ini):
    temp = cells.classProb * cells.ds.area_factor * egamma
    temp = temp.sum(dim='cell_id')
    ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])
    TotPredicted = ClassTotPredicted.drop('Zero', dim='class_name').sum(dim='class_name')
    return TotPredicted


def iss_summary(cells, spots):
    '''
    returns a datafram summarising the main feaures of each cell, ie gene counts and cell types
    :param spots:
    :return:
    '''
    x = cells.ds.x.values
    y = cells.ds.y.values
    cell_id = cells.ds.index.tolist()

    gene_count = cells.geneCount(spots)
    class_prob = cells.classProb
    gene_names = gene_count.gene_name.values
    class_names = class_prob.class_name.values

    tol = 0.001

    logger.info('starting collecting data ...')
    N = len(cell_id)
    isCount_nonZero = [gene_count.values[n, :] > tol for n in range(N)]
    name_list = [gene_names[isCount_nonZero[n]] for n in range(N)]
    count_list = [gene_count[n, isCount_nonZero[n]].values for n in range(N)]

    isProb_nonZero = [class_prob.values[n, :] > tol for n in range(N)]
    class_name_list = [class_names[isProb_nonZero[n]] for n in range(N)]
    prob_list = [class_prob.values[n, isProb_nonZero[n]] for n in range(N)]

    iss_df = pd.DataFrame({'Cell_Num': cells.ds.index.tolist(),
                            'X': cells.ds.x.values,
                            'Y': cells.ds.y.values,
                            'Genenames': name_list,
                            'CellGeneCount': count_list,
                            'ClassName': class_name_list,
                            'Prob': prob_list
                            },
                           columns=['Cell_Num', 'X', 'Y', 'Genenames', 'CellGeneCount', 'ClassName', 'Prob']
                           )
    iss_df.set_index(['Cell_Num'])

    # Sanity check. Only the dummy cell where the misreads are going to be assigned to should not have
    # proper coordinates
    mask = np.isnan(iss_df.X.values) & np.isnan(iss_df.Y.values)
    assert sum(mask) == 1, 'You should have exactly one dummy cell with nan valued coordinates. ' \
                           'It will be used to assign the misreads'

    # Remove that dummy cell from data to be rendered by the viewer
    iss_df = iss_df[~mask]
    logger.info('Data collected!')

    return iss_df


def summary(spots):
    # check for duplicates (ie spots with the same coordinates with or without the same gene name).
    # I dont know why but it can happen. Misfire during spot calling maybe?
    is_duplicate = spots.data.duplicated(subset=['x', 'y'])

    num_rows = spots.data.shape[0]

    cell_prob = spots.call.cell_prob.values
    neighbors = spots.call.neighbors.values
    p = [cell_prob[i, :] for i in range(num_rows)]
    nbrs = [neighbors[i, :] for i in range(num_rows)]
    max_nbrs = [neighbors[i, idx] for i in range(num_rows) for idx in [np.argmax(cell_prob[i, :])]]

    out = pd.DataFrame({'Gene': spots.data.gene_name,
                        'Expt': spots.data.gene_id,
                        'x': spots.data.x,
                        'y': spots.data.y,
                        'neighbour': max_nbrs,
                        'neighbour_array': nbrs,
                        'neighbour_prob': p})

    return out


def collect_data(cells, spots):
    '''
    Collects data for the viewer
    :param cells:
    :param spots:
    :return:
    '''
    iss_df = iss_summary(cells, spots)
    gene_df = summary(spots)
    return iss_df, gene_df
