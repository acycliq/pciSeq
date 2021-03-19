import numpy as np
import pandas as pd
import numpy_groupies as npg
from typing import Tuple
from pciSeq.src.cell_call.datatypes import Cells, Spots, Genes
from pciSeq.src.cell_call.singleCell import sc_expression_data
from pciSeq.src.cell_call.summary import collect_data
import pciSeq.src.cell_call.utils as utils
import gc
import os
import time
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
        self.scData = sc_expression_data(self.genes, scRNAseq, self.config)
        self.nC = self.cells.nC                         # number of cells
        self.nG = self.genes.nG                         # number of genes
        self.nK = len(self.scData.class_name)           # number of classes
        self.nS = self.spots.nS                         # number of spots
        self.nN = self.config['nNeighbors'] + 1         # number of neighbouring cells, candidates for being parent
                                                        # cell of any given spot. The last one will be used for the
                                                        # misread spots.

    def initialise(self):
        # initialise the gene efficiencies and the the starting
        # spot to parent cell assignment
        self.cells.class_names = self.scData.coords['class_name'].values

        self.cells.prior = np.append([.5 * np.ones(self.nK - 1) / self.nK], 0.5)
        self.cells.classProb = np.tile(self.cells.prior, (self.nC, 1))
        self.genes.eta = np.ones(self.nG)
        # self.spots.gamma_bar = np.ones([self.cells.num_cells, self.single_cell_data['mean_expression'].shape[1], self.genes.nG])
        self.cells.alpha = np.ones(self.nC) * self.cells.rho_1 / self.cells.rho_2
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id, self.config)

        # self.cells.centroid_upd(self.spots)
        # self.cells.cov_upd(self.spots)
        # self.spots.init_call(self.cells, self.config)
        self.spots.gamma_bar = np.ones([self.nC, self.nG, self.nK])

        # self.cells.geneCount_upd(self.spots)

    # -------------------------------------------------------------------- #
    def run(self):
        # make a list to keep the main attributes of the fitted bivariate normal (centroids, correlation etc...)
        ellipsoid_attr = [self.cells.ellipsoid_attributes]
        self.initialise()
        # self.spots.init_call(self.cells, self.config)

        p0 = None
        iss_df = None
        gene_df = None
        max_iter = self.config['max_iter']
        for i in range(max_iter):
            ellipsoid_attr.append(self.cells.ellipsoid_attributes)
            # logger.info('initial centroid of cell: 2259 is %s' % self.cells.ini_centroids().iloc[2259].values)
            # logger.info('centroid of cell: 2259 is %s' % self.cells.centroid.iloc[2259].values)
            # logger.info('(mu, rho, sigma_x, sigma_y) of cell: 2259 is %s' % list(self.cells.ellipsoid_attributes[2259]))

            self.geneCount_upd()

            self.alpha_upd()

            # 1. calc expected gamma
            # logger.info('calc expected gamma')
            self.gamma_upd()
            # self.mean_gamma, self.mean_loggamma = self.expected_gamma()

            self.centroid_upd()
            self.cov_upd()

            # 2 assign cells to cell types
            # logger.info('cell to cell type')
            self.prob_cell_to_cellType()

            # 3 assign spots to cells
            # logger.info('spot to cell')
            self.prob_spots_to_cell()

            # 4 update gene efficiency
            # logger.info('update gamma')
            self.eta_upd()


            converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info('Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            if converged:
                iss_df, gene_df = collect_data(self.cells, self.spots, self.genes)
                # print("Success!!")
                break

            if i == max_iter-1:
                logger.info('Loop exhausted. Exiting with convergence status: %s' % converged)
        return iss_df, gene_df

    # -------------------------------------------------------------------- #
    def geneCount_upd(self):
        '''
        Produces a matrix numCells-by-numGenes where element at position (c,g) keeps the expected
        number of gene g  in cell c.
        :param spots:
        :return:
        '''
        start = time.time()

        # nN = self.cells.config['nNeighbors'] + 1
        CellGeneCount = np.zeros([self.nC, self.nG])

        # name = spots.gene_panel.index.values
        spot_id = self.spots.gene_id
        for n in range(self.nN):
            c = self.spots.parent_cell_id[:, n]
            group_idx = np.vstack((c[None, :], spot_id[None, :]))
            a = self.spots.parent_cell_prob[:, n]
            accumarray = npg.aggregate(group_idx, a, func="sum", size=(self.nC, self.nG))
            if n == self.nN - 1:
                # the last neighbouring cell is for the misreads
                self.cells.background_counts = accumarray
            else:
                CellGeneCount = CellGeneCount + accumarray

        end = time.time()
        # print('time in geneCount: ', end - start)
        # CellGeneCount = xr.DataArray(CellGeneCount, coords=[_id, name], dims=['cell_id', 'gene_name'])
        # self.CellGeneCount = CellGeneCount

        # print(self.background_counts.sum())
        # print(CellGeneCount.sum(axis=1).sum())
        # assert self.background_counts.sum() + CellGeneCount.sum(axis=1).sum() == spots.data.shape[0], \
        #     "The sum of the background spots and the cell gene counts should be equal to the total number of spots"
        self.cells.geneCount = CellGeneCount

    # -------------------------------------------------------------------- #
    def gamma_upd(self):
        """
         Implements equation (3) of the Qian paper
        """
        self.cells.alpha
        cells = self.cells
        spots = self.spots
        sc = self.scData
        cfg = self.config
        scaled_mean = np.einsum('c, gk -> cgk', self.cells.alpha, sc.mean_expression)
        # scaled_mean = np.einsum('c, gk -> cgk', cells.cell_props['area_factor'], sc.mean_expression)
        rho = cfg['rSpot'] + cells.geneCount
        beta = cfg['rSpot'] + scaled_mean

        expected_gamma = utils.gammaExpectation(rho, beta)
        expected_loggamma = utils.logGammaExpectation(rho, beta)

        del rho
        del beta
        gc.collect()
        del gc.garbage[:]
        self.spots.gamma_bar = expected_gamma
        self.spots.log_gamma_bar = expected_loggamma
        # return expected_gamma, expected_loggamma

    # -------------------------------------------------------------------- #
    def prob_cell_to_cellType(self):
        """
        return a an array of size numCells-by-numCellTypes where element in position [i,j]
        keeps the probability that cell i has cell type j
        Implements equation (2) of the Qian paper
        :param spots:
        :param config:
        :return:
        """

        # gene_gamma = self.genes.eta
        sc = self.scData
        ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.alpha, self.genes.eta, sc.mean_expression) + self.config['SpotReg']
        # ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.cell_props['area_factor'], self.genes.eta, sc.mean_expression) + self.config['SpotReg']
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cells.log_prior
        pCellClass = utils.softmax(wCellClass, axis=1)

        ## self.cells.classProb = pCellClass
        # logger.info('Cell 0 is classified as %s with prob %4.8f' % (
        #     self.cells.class_names[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
        # logger.info('cell ---> cell class probabilities updated')
        self.cells.classProb = pCellClass
        # return pCellClass

    # -------------------------------------------------------------------- #
    def prob_spots_to_cell(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        nN = self.config['nNeighbors'] + 1
        nS = self.spots.data.gene_name.shape[0]

        # initialise array with the misread density
        aSpotCell = np.zeros([nS, nN]) + np.log(self.config['MisreadDensity'])
        gn = self.spots.data.gene_name.values
        expected_spot_counts = self.scData['log_mean_expression'].sel({'gene_name': gn}).data

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        for n in range(nN - 1):
            spots_name = self.spots.data.gene_name.values
            self.scData['log_mean_expression'].sel({'gene_name': spots_name})

            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            # cp = self.cells.classProb.sel({'cell_id': sn}).data
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            # term_1 = (expected_spot_counts * cp).sum(axis=1)
            term_1 = np.einsum('ij, ij -> i', expected_spot_counts, cp)

            # logger.info('genes.spotNo should be something like spots.geneNo instead!!')
            expectedLog = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]
            # expectedLog = utils.bi2(self.elgamma.data, [nS, nK], sn[:, None], self.spots.data.gene_id.values[:, None])

            term_2 = np.einsum('ij, ij -> i', cp, expectedLog)  # same as np.sum(cp * expectedLog, axis=1) but bit faster

            loglik = self.spots.mvn_loglik(self.spots.xy_coords, sn, self.cells)
            aSpotCell[:, n] = term_1 + term_2 + loglik
            # logger.info('')
        wSpotCell = aSpotCell # + self.spots.loglik(self.cells, self.config)

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.parent_cell_prob = pSpotNeighb

        self.geneCount_upd()
        # self.spots.update_cell_prob(pSpotNeighb, self.cells)
        # logger.info('spot ---> cell probabilities updated')

    # -------------------------------------------------------------------- #
    def eta_upd(self):
        """
        Calcs the expected eta
        Implements equation (5) of the Qian paper
        """
        grand_total = self.cells.background_counts.sum() + self.cells.total_counts.sum()
        assert round(grand_total) == self.spots.data.shape[0], \
            'The sum of the background spots and the total gene counts should be equal to the number of spots'

        zero_prob = self.cells.classProb[:, -1] # probability a cell being a zero expressing cell
        zero_class_counts = self.spots.zero_class_counts(self.spots.gene_id, zero_prob)
        class_total_counts = self.cells.geneCountsPerKlass(self.scData, self.spots.gamma_bar, self.config)

        # TotPredictedB = np.bincount(self.spots.gene_id, self.spots.adj_cell_prob[:, -1].data)
        # assert np.all(TotPredictedB == self.cells.background_counts)
        background_counts = self.cells.background_counts
        nom = self.config['rGene'] + self.spots.counts_per_gene - background_counts - zero_class_counts
        denom = self.config['rGene'] + class_total_counts
        res = nom / denom

        # Finally, update gene_gamma
        self.genes.eta = res.values

    # -------------------------------------------------------------------- #
    def alpha_upd(self):
        ## manual way to do the denom, useful for debugging
        ## Pick up a cell, lets say the first one and get the ids of its most likely cell types
        # mask = np.argsort(-1 * self.classProb[1, :])

        ## select the top 3 most likely cell types
        # mask = mask[: 3]
        ## find the theoretical gene counts that this particular cell should have
        # out = np.einsum('cg, g, gc, c -> ', gamma_bar[1, :, mask], genes.eta, sc.mean_expression.data[:, mask], self.classProb[1, mask])

        N_c = self.cells.total_counts
        zeta_bar = self.cells.classProb
        mu = self.scData['mean_expression'].data   # need to add constant too!

        denom = np.einsum('gk, ck, cgk, g -> c', mu, zeta_bar, self.spots.gamma_bar, self.genes.eta)
        alpha_bar = (N_c + self.cells.rho_1) / (denom + self.cells.rho_2)
        logger.info('alpha range: %s ' % [alpha_bar.min(), alpha_bar.max()])
        self.cells.alpha = alpha_bar

    # -------------------------------------------------------------------- #
    def centroid_upd(self):
        spots = self.spots

        # get the total gene counts per cell
        N_c = self.cells.total_counts

        xy_spots = spots.xy_coords
        prob = spots.parent_cell_prob
        n = self.cells.config['nNeighbors'] + 1

        # mulitply the x coord of the spots by the cell prob
        a = np.tile(xy_spots[:, 0], (n, 1)).T * prob

        # mulitply the y coord of the spots by the cell prob
        b = np.tile(xy_spots[:, 1], (n, 1)).T * prob

        # aggregated x and y coordinate
        idx = spots.parent_cell_id
        x_agg = npg.aggregate(idx.ravel(), a.ravel(), size=len(N_c))
        y_agg = npg.aggregate(idx.ravel(), b.ravel(), size=len(N_c))

        # x_agg = np.zeros(N_c.shape)
        # mask = np.arange(len(_x_agg))
        # x_agg[mask] = _x_agg
        #
        # y_agg = np.zeros(N_c.shape)
        # mask = np.arange(len(_y_agg))
        # y_agg[mask] = _y_agg

        # get the estimated cell centers
        x_bar = np.nan * np.ones(N_c.shape)
        y_bar = np.nan * np.ones(N_c.shape)
        x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
        y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]
        # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
        xy_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T)))

        # if you have a value for the estimated centroid use that, otherwise
        # use the initial (starting values) centroids
        ini_cent = self.cells.ini_centroids()
        xy_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y']])))

        # # sanity check. NaNs or Infs should appear together
        # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
        # use the fitted centroids where possible otherwise use the initial ones
        xy_bar[np.isfinite(x_bar)] = xy_bar_fitted[np.isfinite(x_bar)]
        self.cells.centroid = pd.DataFrame(xy_bar, columns=['x', 'y'])
        # print(np.array(list(zip(x_bar.T, y_bar.T))))

    # -------------------------------------------------------------------- #
    def cov_upd(self):
        spots = self.spots

        # first get the scatter matrix
        S = self.cells.scatter_matrix(spots)  # sample sum of squares
        cov_0 = self.cells.ini_cov()
        nu_0 = self.cells.nu_0
        S_0 = cov_0 * nu_0 # prior sum of squarea
        N_c = self.cells.total_counts
        d = 2
        denom = np.ones([self.nC, 1, 1])
        denom[:, 0, 0] = N_c + nu_0
        # Note: need to add code to handle the case N_c + nu_0 <= d + 2
        cov = (S + S_0) / denom

        # delta = self.cells.ledoit_wolf(self.spots, cov)
        # stein = self.cells.stein(cov)

        delta = 0.5
        # logger.info('Mean shrinkage %4.2f' % delta.mean())
        # logger.info('cell 601 shrinkage %4.2f, %4.2f' % (shrinkage[601], sh[601]))
        logger.info('cell 601 gene counts %d' % self.cells.total_counts[601])
        # logger.info('cell 605 shrinkage %4.2f, %4.2f' % (shrinkage[605], sh[605]))
        logger.info('cell 605 gene counts %d' % self.cells.total_counts[605])
        # logger.info('cell 610 shrinkage %4.2f, %4.2f' % (shrinkage[610], sh[610]))
        logger.info('cell 610 gene counts %d' % self.cells.total_counts[610])

        # delta = delta.reshape(self.nC, 1, 1)
        cov = delta*cov_0 + (1-delta)*cov
        self.cells.cov = cov






