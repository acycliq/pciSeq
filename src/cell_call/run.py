from src.cell_call.systemData import Cells, Spots, Prior
from src.cell_call.singleCell import sc_expression_data
import src.cell_call.common as common
import src.cell_call.utils as utils
import os
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


class VarBayes:
    def __init__(self, config):
        self.config = config
        self.cells = Cells(self.config)
        self.spots = Spots(self.config)
        self.single_cell_data = sc_expression_data(self.spots, self.config)

    def run(self):
        class_names = self.single_cell_data.coords['class_name'].values
        prior = Prior(class_names)
        self.spots.init_call(self.cells, self.config)

        p0 = None
        iss_df = None
        gene_df = None
        for i in range(self.config['max_iter']):
            # 1. calc expected gamma
            logger.info('calc expected gamma')
            egamma, elgamma = common.expected_gamma(self.cells, self.spots, self.single_cell_data, self.config)

            # 2 assign cells to cell types
            logger.info('cell to cell type')
            common.celltype_assignment(self.cells, self.spots, prior, self.single_cell_data, self.config)

            # 3 assign spots to cells
            logger.info('spot to cell')
            common.call_spots(self.spots, self.cells, self.single_cell_data, prior, elgamma, self.config)

            # 4 update eta
            logger.info('update gamma')
            common.updateGamma(self.cells, self.spots, self.single_cell_data, egamma, self.config)

            converged, delta = utils.isConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info('Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.call.cell_prob.values

            if converged:
                # cells.iss_summary(spots)
                # spots.summary()
                iss_df, gene_df = common.collect_data(self.cells, self.spots)
                print("Success!!")
                break

        return iss_df, gene_df




