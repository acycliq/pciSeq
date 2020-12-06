from src.cell_call.systemData import Cells, Spots, Prior
from src.cell_call.singleCell import sc_expression_data
import src.cell_call.common as common
import pandas as pd
import src.cell_call.utils as utils
import os
import logging


dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


def varBayes(my_config):
    saFile = os.path.join(dir_path, my_config['saFile'])
    cells = Cells(my_config)

    logger.info('********* Getting spot attributes from %s **********', saFile)
    sa_df = pd.read_csv(saFile)
    if my_config['drop_nan']:
        sa_df = sa_df.dropna()  ##  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE  CAREFUL HERE
        logger.info('**** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs ***** I HAVE REMOVED NaNs****')

    sa_df = sa_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

    # remove a gene if it is on the exclude list
    gene_mask = [True if d not in my_config['exclude_genes'] else False for d in sa_df.target]
    sa_df = sa_df.loc[gene_mask]

    spots = Spots(sa_df)
    single_cell_data = sc_expression_data(spots, my_config)
    prior = Prior(single_cell_data.coords['class_name'].values)
    spots.init_call(cells, my_config)

    p0 = None
    iss_df = None
    gene_df = None
    for i in range(my_config['max_iter']):
        # 1. calc expected gamma
        logger.info('calc expected gamma')
        egamma, elgamma = common.expected_gamma(cells, spots, single_cell_data, my_config)

        # 2 call cells
        logger.info('cell to cell type')
        common.celltype_assignment(cells, spots, prior, single_cell_data, my_config)

        # 3 call spots
        logger.info('spot to cell')
        common.call_spots(spots, cells, single_cell_data, prior, elgamma, my_config)

        # 4 update eta
        logger.info('update gamma')
        common.updateGamma(cells, spots, single_cell_data, egamma, my_config)

        converged, delta = utils.isConverged(spots, p0, my_config['CellCallTolerance'])
        logger.info('Iteration %d, mean prob change %f' % (i, delta))

        # replace p0 with the latest probabilities
        p0 = spots.call.cell_prob.values

        if converged:
            # cells.iss_summary(spots)
            # spots.summary()
            iss_df, gene_df = common.collect_data(cells, spots)
            print("Success!!")
            break

    return iss_df, gene_df


logger.info('Done')

