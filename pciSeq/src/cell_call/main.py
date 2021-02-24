import numpy as np
from pciSeq.src.cell_call.datatypes import Cells, Spots, Genes
from pciSeq.src.cell_call.singleCell import sc_expression_data
from pciSeq.src.cell_call.summary import collect_data
import pciSeq.src.cell_call.utils as utils
import gc
import os
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()

# Note: 21-Feb-2020
# I think there are TWO  bugs.
#
# ********** BUG 1 **********
# I should be using: mean_expression = mean_expression + config['SpotReg']
# instead of adding the SpotReg value every time I am calling mean_expression.
# The ScaledExp calculation, in the code below, for example is not correct because the config['SpotReg'] constant
# should not be added at the very end but should be moved inside the parentheses and added directly to
# the sc.mean_expression array
#
# ********** BUG 2 **********
# Gene efficienty (eta) is not also handled properly. For eta to have a gamma distribution gamma(r, r/eta_0) with
# mean=0.2 and variance 20 the hyperparameter should be r = 0.002 and r/eta_0 = 0.002/0.2 = 0.01.
# The mean value (ie 0.2) is applied on the single cell gene expression right after we have calculated the mean
# class expression.
# Then, in the rest of the code, gene efficiency is initialised as a vector of ones. While the vector eta is
# getting multiplied with the mean expression counts when we derive the cell to cell class assignment, we do not
# do the same when we calculate the gamma (equation 3 in the NMETH paper. We multiply the mean class expression
# by the cell area factor and the gene efficiency is missing (see ScaledMean in matlab code)
# I think It would better to:
# NOT multiply the mean class expressions counts by 0.2 and the initialise eta as a vector of ones
# BUT instead:
# do not scale down the single cell gene expression data by 0.2 at the very beginning
# initialise eta as a vector 0.2*np.ones([Ng, 1])
# include eta in the calculations equations 2 and 3 and

class VarBayes:
    def __init__(self, _cells_df, _spots_df, scRNAseq, config):
        self.config = config
        self.cells = Cells(_cells_df, config)
        self.spots = Spots(_spots_df, config)
        self.genes = Genes(self.spots)
        self.single_cell_data = sc_expression_data(self.genes, scRNAseq, self.config)
        class_names = self.single_cell_data.coords['class_name'].values
        self.cells.log_prior = self.cells.prior(class_names)
        self.cells.class_names = class_names

        self.egamma = None
        self.elgamma = None

    def initialise(self):
        # initialise the gene efficiencies and the the starting
        # spot to parent cell assignment
        # No need to initialise the spot gammas, they are the very first
        # arrays to be calculated inside the run () method
        self.genes.eta = np.ones(self.genes.nG)
        self.spots.init_call(self.cells, self.config)
        # self.cells.geneCount_upd(self.spots)

    # -------------------------------------------------------------------- #
    def run(self):
        self.initialise()
        # self.spots.init_call(self.cells, self.config)

        p0 = None
        iss_df = None
        gene_df = None
        max_iter = self.config['max_iter']
        for i in range(max_iter):
            # 1. calc expected gamma
            # logger.info('calc expected gamma')
            self.mean_gamma, self.mean_loggamma = self.expected_gamma()

            # 2 assign cells to cell types
            # logger.info('cell to cell type')
            self.cells.classProb = self.update_classProb()

            # 3 assign spots to cells
            # logger.info('spot to cell')
            self.assign_spots()

            # 4 update gene efficiency
            # logger.info('update gamma')
            self.update_eta()

            converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info('Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.adj_cell_prob

            if converged:
                iss_df, gene_df = collect_data(self.cells, self.spots, self.genes)
                # print("Success!!")
                break

            if i == max_iter-1:
                logger.info('Loop exhausted. Exiting with convergence status: %s' % converged)
        return iss_df, gene_df

    # -------------------------------------------------------------------- #
    def expected_gamma(self):
        cells = self.cells
        spots = self.spots
        sc = self.single_cell_data
        cfg = self.config
        scaled_mean = np.einsum('c, gk -> cgk', cells.cell_props['area_factor'], sc.mean_expression)
        rho = cfg['rSpot'] + cells.geneCount
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

        # gene_gamma = self.genes.eta
        sc = self.single_cell_data
        ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.cell_props['area_factor'], self.genes.eta, sc.mean_expression) + self.config['SpotReg']
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cells.log_prior
        pCellClass = utils.softmax(wCellClass, axis=1)

        ## self.cells.classProb = pCellClass
        # logger.info('Cell 0 is classified as %s with prob %4.8f' % (
        #     self.cells.class_names[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
        # logger.info('cell ---> cell class probabilities updated')
        return pCellClass

    # -------------------------------------------------------------------- #
    def assign_spots(self):
        # spot to cell assignment
        nN = self.config['nNeighbors'] + 1
        nS = self.spots.data.gene_name.shape[0]
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
            expectedLog = self.mean_loggamma[self.spots.adj_cell_id[:, n], self.spots.gene_id]
            # expectedLog = utils.bi2(self.elgamma.data, [nS, nK], sn[:, None], self.spots.data.gene_id.values[:, None])

            term_2 = np.einsum('ij,ij -> i', cp, expectedLog)  # same as np.sum(cp * expectedLog, axis=1) but bit faster
            aSpotCell[:, n] = term_1 + term_2
        wSpotCell = aSpotCell + self.spots.loglik(self.cells, self.config)

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.adj_cell_prob = pSpotNeighb, self.cells
        # logger.info('spot ---> cell probabilities updated')

    # -------------------------------------------------------------------- #
    def update_eta(self):
        grand_total = self.cells.background_counts.sum() + self.cells.total_counts.sum()
        assert round(grand_total) == self.spots.data.shape[0], \
            'The sum of the background spots and the total gene counts should be equal to the number of spots'

        zero_prob = self.cells.classProb[:, -1] # probability a cell being a zero expressing cell
        zero_class_counts = self.spots.zero_class_counts(self.spots.gene_id, zero_prob)
        class_total_counts = self.cells.geneCountsPerKlass(self.single_cell_data, self.mean_gamma, self.config)

        # TotPredictedB = np.bincount(self.spots.gene_id, self.spots.adj_cell_prob[:, -1].data)
        # assert np.all(TotPredictedB == self.cells.background_counts)
        background_counts = self.cells.background_counts
        nom = self.config['rGene'] + self.spots.counts_per_gene - background_counts - zero_class_counts
        denom = self.config['rGene'] + class_total_counts
        res = nom / denom

        # Finally, update gene_gamma
        self.genes.eta = res.values







