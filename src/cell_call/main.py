import numpy as np
from src.cell_call.datatypes import Cells, Spots, Genes
from src.cell_call.singleCell import sc_expression_data
from src.cell_call.summary import collect_data
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
        self.genes = Genes(self.spots)
        self.single_cell_data = sc_expression_data(self.genes, self.config)
        class_names = self.single_cell_data.coords['class_name'].values
        self.cells.log_prior = self.cells.prior(class_names)
        self.cells.class_names = class_names

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
            self.assign_spots()

            # 4 update gene efficiency
            logger.info('update gamma')
            self.update_eta()

            converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info('Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.adj_cell_prob

            if converged:
                iss_df, gene_df = collect_data(self.cells, self.spots, self.genes)
                print("Success!!")
                break
        return iss_df, gene_df

    # -------------------------------------------------------------------- #
    def expected_gamma(self):
        cells = self.cells
        spots = self.spots
        sc = self.single_cell_data
        cfg = self.config
        scaled_mean = np.einsum('c, gk -> cgk', cells.cell_props['area_factor'], sc.mean_expression)
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
        """
        return a an array of size numCells-by-numCellTypes where element in position [i,j]
        keeps the probability that cell i has cell type j
        :param spots:
        :param config:
        :return:
        """

        gene_gamma = self.genes.gamma
        sc = self.single_cell_data
        ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.cell_props['area_factor'], gene_gamma, sc.mean_expression) + self.config['SpotReg']
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount(self.spots)
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cells.log_prior
        pCellClass = utils.softmax(wCellClass, axis=1)

        # self.cells.classProb = pCellClass
        logger.info('Cell 0 is classified as %s with prob %4.8f' % (
            self.cells.class_names[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
        logger.info('cell ---> cell class probabilities updated')
        return pCellClass

    # -------------------------------------------------------------------- #
    def assign_spots(self):
        # spot to cell assignment
        nN = self.config['nNeighbors'] + 1
        nS = self.spots.data.gene_name.shape[0]
        # nK = self.cell_prior.nK
        aSpotCell = np.zeros([nS, nN])
        gn = self.spots.data.gene_name.values
        expected_spot_counts = self.single_cell_data['log_mean_expression'].sel({'gene_name': gn}).data
        for n in range(nN - 1):
            spots_name = self.spots.data.gene_name.values
            self.single_cell_data['log_mean_expression'].sel({'gene_name': spots_name})

            # get the spots' nth-closest cell
            sn = self.spots.adj_cell_id[:, n]

            # get the respective cell type probabilities
            # cp = self.cells.classProb.sel({'cell_id': sn}).data
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            # term_1 = (expected_spot_counts * cp).sum(axis=1)
            term_1 = np.einsum('ij, ij -> i', expected_spot_counts, cp)

            # logger.info('genes.spotNo should be something like spots.geneNo instead!!')
            expectedLog = self.elgamma[self.spots.adj_cell_id[:, n], self.spots.gene_id]
            # expectedLog = utils.bi2(self.elgamma.data, [nS, nK], sn[:, None], self.spots.data.gene_id.values[:, None])

            term_2 = np.einsum('ij,ij -> i', cp, expectedLog)  # same as np.sum(cp * expectedLog, axis=1) but bit faster
            aSpotCell[:, n] = term_1 + term_2
        wSpotCell = aSpotCell + self.spots.loglik(self.cells, self.config)

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.adj_cell_prob = pSpotNeighb
        logger.info('spot ---> cell probabilities updated')

    # -------------------------------------------------------------------- #
    def update_eta(self):
        # Maybe I should rename that to eta (not gamma). In the paper the symbol is eta
        zero_prob = self.cells.classProb[:, -1] # probability a cell being a zero expressing cell
        TotPredictedZ = self.spots.TotPredictedZ(self.spots.gene_id, zero_prob)

        TotPredicted = self.cells.geneCountsPerKlass(self.single_cell_data, self.egamma, self.config)

        TotPredictedB = np.bincount(self.spots.gene_id, self.spots.adj_cell_prob[:, -1].data)

        nom = self.config['rGene'] + self.spots.counts_per_gene - TotPredictedB - TotPredictedZ
        denom = self.config['rGene'] + TotPredicted
        res = nom / denom

        # Finally, update gene_gamma
        self.genes.gamma = res.values







