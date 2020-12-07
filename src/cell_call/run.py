import numpy as np
import xarray as xr
from src.cell_call.systemData import Cells, Spots, Cell_prior
from src.cell_call.singleCell import sc_expression_data
import src.cell_call.common as common
import src.cell_call.utils as utils
import gc
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
        class_names = self.single_cell_data.coords['class_name'].values
        self.cell_prior = Cell_prior(class_names)

        self.egamma = None
        self.elgamma = None

    # -------------------------------------------------------------------- #
    def run(self):
        self.spots.init_call(self.cells, self.config)

        p0 = None
        iss_df = None
        gene_df = None
        for i in range(self.config['max_iter']):
            # 1. calc expected gamma
            logger.info('calc expected gamma')
            self.egamma, self.elgamma = self.expected_gamma()

            # 2 assign cells to cell types
            logger.info('cell to cell type')
            self.cells.classProb = self.update_classProb()

            # 3 assign spots to cells
            logger.info('spot to cell')
            self.call_spots()

            # 4 update eta
            logger.info('update gamma')
            self.update_eta()

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

    # -------------------------------------------------------------------- #
    def expected_gamma(self):
        cells = self.cells
        spots = self.spots
        sc = self.single_cell_data
        cfg = self.config
        scaled_mean = cells.cell_props.area_factor.to_xarray() * sc.mean_expression
        rho = cfg['rSpot'] + cells.geneCount(spots)
        beta = cfg['rSpot'] + scaled_mean

        expected_gamma = utils.gammaExpectation(rho, beta)
        expected_loggamma = utils.logGammaExpectation(rho, beta)

        del rho
        del beta
        gc.collect()
        del gc.garbage[:]

        return expected_gamma, expected_loggamma

    # -------------------------------------------------------------------- #
    def update_classProb(self):
        '''
        return a an array of size numCells-by-numCellTypes where element in position [i,j]
        keeps the probability that cell i has cell type j
        :param spots:
        :param config:
        :return:
        '''

        gene_gamma = self.spots.gene_panel.gene_gamma
        sc = self.single_cell_data
        ScaledExp = self.cells.cell_props.area_factor.to_xarray() * gene_gamma.to_xarray() * sc.mean_expression + self.config['SpotReg']
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount(self.spots)
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cell_prior.logvalue
        pCellClass = xr.apply_ufunc(utils.softmax_nolan, wCellClass, 1, 1)

        # self.cells.classProb = pCellClass
        logger.info('Cell 0 is classified as %s with prob %4.8f' % (
            self.cell_prior.name[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
        logger.info('cell ---> klass probabilities updated')
        return pCellClass

    # -------------------------------------------------------------------- #
    def call_spots(self):
        # spot to cell assignment
        nN = self.spots.call.neighbors.shape[1]
        nS = self.spots.data.gene_name.shape[0]
        nK = self.cell_prior.nK
        aSpotCell = np.zeros([nS, nN])
        gn = self.spots.data.gene_name.values
        expected_spot_counts = self.single_cell_data['log_mean_expression'].sel({'gene_name': gn}).data
        for n in range(nN - 1):
            spots_name = self.spots.data.gene_name.values
            self.single_cell_data['log_mean_expression'].sel({'gene_name': spots_name})

            # get the spots' nth-closest cell
            sn = self.spots.call.neighbors.loc[:, n].values

            # get the respective cell type probabilities
            cp = self.cells.classProb.sel({'cell_id': sn}).data

            # multiply and sum over cells
            term_1 = (expected_spot_counts * cp).sum(axis=1)

            # logger.info('genes.spotNo should be something line spots.geneNo instead!!')
            # expectedLog = elgamma.data[spots.call.neighbors[:, n], spots.data.gene_id]
            expectedLog = utils.bi2(self.elgamma.data, [nS, nK], sn[:, None], self.spots.data.gene_id.values[:, None])
            term_2 = np.sum(cp * expectedLog, axis=1)
            aSpotCell[:, n] = term_1 + term_2
        wSpotCell = aSpotCell + self.spots.loglik(self.cells, self.config)

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax2(wSpotCell)
        self.spots.call.cell_prob.data = pSpotNeighb
        logger.info('spot ---> cell probabilities updated')

    # -------------------------------------------------------------------- #
    def update_eta(self):
        # Maybe I should rename that to eta (not gamma). In the paper the symbol is eta
        TotPredictedZ = self.spots.TotPredictedZ(self.spots.data.gene_id.values, self.cells.classProb.sel({'class_name': 'Zero'}).data)

        TotPredicted = self.cells.geneCountsPerKlass(self.single_cell_data, self.egamma, self.config)

        TotPredictedB = np.bincount(self.spots.data.gene_id.values, self.spots.call.cell_prob.loc[:, 3].data)

        nom = self.config['rGene'] + self.spots.gene_panel.total_spots - TotPredictedB - TotPredictedZ
        denom = self.config['rGene'] + TotPredicted
        res = nom / denom

        # Finally, update gene_gamma
        self.spots.gene_panel.gene_gamma[res.index] = res.values
        # cells.expectedGamma = nom / denom






