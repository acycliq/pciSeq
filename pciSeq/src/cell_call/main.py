import numpy as np
import pandas as pd
import numpy_groupies as npg
import datetime
import time
import scipy.spatial as spatial
from pciSeq.src.cell_call.datatypes import Cells, Spots, Genes, SingleCell, CellType
from pciSeq.src.cell_call.summary import collect_data
import pciSeq.src.cell_call.utils as utils
from pciSeq.src.diagnostics.utils import redis_db
from pciSeq.src.diagnostics.launch_diagnostics import launch_dashboard
# from pciSeq.src.diagnostics.launch_diagnostics_dummy import launch_dashboard
import pciSeq.src.diagnostics.config as diagnostics_cfg
from multiprocessing.dummy import Pool as ThreadPool
from pciSeq.src.cell_call.log_config import logger


class VarBayes(object):
    def __init__(self, _cells_df, _spots_df, scRNAseq, config):
        self.config = config
        self.redis_db = redis_db(flush=True)
        self.cells = Cells(_cells_df, config)
        self.spots = Spots(_spots_df, config)
        self.genes = Genes(self.spots)
        self.single_cell = SingleCell(scRNAseq, self.genes.gene_panel, self.config)
        self.cellTypes = CellType(self.single_cell, config)
        self.cells.class_names = self.single_cell.classes
        self.nC = self.cells.nC  # number of cells
        self.nG = self.genes.nG  # number of genes
        self.nK = self.cellTypes.nK  # number of classes
        self.nS = self.spots.nS  # number of spots
        self.nN = self.config['nNeighbors'] + 1     # number of the closest nearby cells, candidates for being parent
                                                    # cell of any given spot. The last cell will be used for the
                                                    # misread spots. (ie cell at position nN is the background)
        self.iter_num = None
        self.has_converged = False

    def initialise(self):
        self.cellTypes.ini_prior()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.genes.init_eta(1, 1 / self.config['Inefficiency'])
        self.spots.parent_cell_id, _ = self.spots.cells_nearby(self.cells)
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id,  # assign a spot to a cell if
                                                              self.config)                # it is within its boundaries
        self.spots.gamma_bar = np.ones([self.nC, self.nG, self.nK]).astype(self.config['dtype'])
        if self.config['launch_diagnostics']:
            logger.info('Launching the diagnostics dashboard')
            launch_dashboard()

    def run(self):
        p0 = None
        iss_df = None
        gene_df = None
        max_iter = self.config['max_iter']

        self.initialise()
        for i in range(max_iter):
            self.iter_num = i

            # 1. For each cell, calc the expected gene counts
            logger.info('1. geneCount_upd')

            self.geneCount_upd()

            # 2. calc expected gamma
            logger.info('2. gamma_upd')
            self.gamma_upd()

            if self.config['relax_segmentation'] or self.config['is_3D']:
                logger.info('3. gaussian_upd')
                self.gaussian_upd()

            # 3. assign cells to cell types
            logger.info('4. cell_to_cellType')
            self.cell_to_cellType()

            # update cell type weights
            # if self.config['relax_segmentation'] or self.config['is_3D']:
            #     self.dalpha_upd()

            # 4. assign spots to cells
            logger.info('6. spots_to_cell')
            # self.spots_to_cell()
            self.spots_to_cell_par()

            # 5. update gene efficiency
            logger.info('7. eta_upd')
            self.eta_upd()

            # 6. Update single cell data
            if self.single_cell.isMissing:
                logger.info('8. mu_upd')
                self.mu_upd()

            self.has_converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info(' Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            logger.info('start db save')
            # self.db_save()
            self.redis_save()
            logger.info('end db save')
            if self.has_converged:
                # self.db_save()
                self.counts_within_radius(20)
                iss_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell, self.config)
                break

            if i == max_iter - 1:
                logger.info(' Loop exhausted. Exiting with convergence status: %s' % self.has_converged)

        return iss_df, gene_df

    # -------------------------------------------------------------------- #
    def geneCount_upd(self):
        """
        Produces a matrix numCells-by-numGenes where element at position (c,g) keeps the expected
        counts of gene g  in cell c.
        """
        # make an array nS-by-nN and fill it with the spots id
        gene_ids = np.tile(self.spots.gene_id, (self.nN, 1)).T

        # flatten it
        gene_ids = gene_ids.ravel()

        # make corresponding arrays for cell_id and probs
        cell_ids = self.spots.parent_cell_id.ravel()
        probs = self.spots.parent_cell_prob.ravel()

        # make the array to be used as index in the group-by operation
        group_idx = np.vstack((cell_ids, gene_ids))

        # For each cell aggregate the number of spots from the same gene.
        # It will produce an array of size nC-by-nG where the entry at (c,g)
        # is the gene counts of gene g within cell c
        N_cg = npg.aggregate(group_idx, probs, size=(self.nC, self.nG))

        # assert N_cg.sum() == self.spots.data.shape[0], \
        #     "The sum of the background spots and the cell gene counts should be equal to the total number of spots"

        # make output. This part needs to be rewritten
        out = np.zeros([self.nC, self.nG])
        out[1:, :] = N_cg[1:, :]

        # cell at position zero is the background
        self.cells.background_counts = N_cg[0, :]
        # Actual cells are on non-zero positions
        self.cells.geneCount = out

    # -------------------------------------------------------------------- #
    def gamma_upd(self):
        """
         Implements equation (3) of the Qian paper
        """
        cells = self.cells
        cfg = self.config
        dtype = self.config['dtype']
        beta = np.einsum('c, gk, g -> cgk',
                         cells.ini_cell_props['area_factor'],
                         self.single_cell.mean_expression,
                         self.genes.eta_bar).astype(dtype) + cfg['rSpot']
        rho = cfg['rSpot'] + cells.geneCount

        self.spots.gamma_bar = self.spots.gammaExpectation(rho, beta)
        self.spots.log_gamma_bar = self.spots.logGammaExpectation(rho, beta)

    # -------------------------------------------------------------------- #
    def cell_to_cellType(self):
        """
        return an array of size numCells-by-numCellTypes where element in position [i,j]
        keeps the probability that cell i has cell type j
        Implements equation (2) of the Qian paper
        :param spots:
        :param config:
        :return:
        """

        # gene_gamma = self.genes.eta
        dtype = self.config['dtype']
        # ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.alpha, self.genes.eta, sc.mean_expression.data) + self.config['SpotReg']
        ScaledExp = np.einsum('c, g, gk -> cgk',
                              self.cells.ini_cell_props['area_factor'],
                              self.genes.eta_bar,
                              self.single_cell.mean_expression
                              ).astype(dtype)
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cellTypes.log_prior
        # if self.config['is_3D']:
        #     wCellClass = np.sum(contr, axis=1) + self.cellTypes.logpi_bar
        # else:
        #     wCellClass = np.sum(contr, axis=1) + np.log(self.cellTypes.prior)
        pCellClass = utils.softmax(wCellClass, axis=1)

        ## self.cells.classProb = pCellClass
        # logger.info('Cell 0 is classified as %s with prob %4.8f' % (
        #     self.cells.class_names[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
        # logger.info('cell ---> cell class probabilities updated')
        self.cells.classProb = pCellClass
        # return pCellClass

    # -------------------------------------------------------------------- #
    # def spots_to_cell_2D_XXX(self):
    #     """
    #     spot to cell assignment.
    #     Implements equation (4) of the Qian paper
    #     """
    #     nN = self.nN
    #     nS = self.spots.data.gene_name.shape[0]
    #
    #     # initialise array with the misread density
    #     aSpotCell = np.zeros([nS, nN])
    #     gn = self.spots.data.gene_name.values
    #     expected_counts = self.single_cell.log_mean_expression.loc[gn].values
    #
    #     # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
    #     my_D = np.zeros([nS, nN])
    #     for n in range(nN - 1):
    #         # get the spots' nth-closest cell
    #         sn = self.spots.parent_cell_id[:, n]
    #
    #         # get the respective cell type probabilities
    #         cp = self.cells.classProb[sn]
    #
    #         # multiply and sum over cells
    #         term_1 = np.einsum('ij, ij -> i', expected_counts, cp)
    #
    #         # logger.info('genes.spotNo should be something like spots.geneNo instead!!')
    #         expectedLog = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]
    #
    #         term_2 = np.einsum('ij, ij -> i', cp, expectedLog)
    #         aSpotCell[:, n] = term_1 + term_2
    #         # my_D[:, n] = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells)
    #         my_covs = self.cells.cov[sn] * np.diag([1, 1, 1])
    #
    #         # my_D[:, n] = self.spots.multiple_logpdfs(self.spots.xyz_coords[:, :2],  self.cells.centroid.values[sn][:, :2], my_covs[:,:2,:2])
    #         my_D[:, n] = self.spots.multiple_logpdfs(self.spots.xyz_coords, self.cells.centroid.values[sn], my_covs)
    #     wSpotCell = aSpotCell + self.spots.loglik(self.cells, self.config)
    #
    #     # update the prob a spot belongs to a neighboring cell
    #     pSpotNeighb = utils.softmax(wSpotCell, axis=1)
    #     self.spots.parent_cell_prob = pSpotNeighb
    #     # Since the spot-to-cell assignments changed you need to update the gene counts now
    #     self.geneCount_upd()
    #     # self.spots.update_cell_prob(pSpotNeighb, self.cells)
    #     # logger.info('spot ---> cell probabilities updated')

    # -------------------------------------------------------------------- #
    def spots_to_cell(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        # nN = self.config['nNeighbors'] + 1
        nN = self.nN
        nS = self.spots.data.gene_name.shape[0]

        misread_density_adj_factor = -1 * np.log(2 * np.pi * self.cells.mcr ** 2) / 2
        misread_density_adj = np.exp(misread_density_adj_factor) * self.config['MisreadDensity']
        # initialise array with the misread density
        wSpotCell = np.zeros([nS, nN]) + np.log(misread_density_adj)
        gn = self.spots.data.gene_name.values
        sc_mean = self.single_cell.log_mean_expression.loc[gn].values
        logeta_bar = self.genes.logeta_bar[self.spots.gene_id]

        logger.info('start loop')
        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        my_D = np.zeros([nS, nN])
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            term_1 = np.einsum('ij, ij -> i', sc_mean, cp)

            log_gamma_bar = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)
            # logger.info('start mvn_loglik')
            loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells)
            # logger.info('end mvn_loglik')
            my_D[:, n] = loglik

            # if 2D, add the bonus for the spots inside the boundaries.
            if not self.config['is_3D']:
                loglik = loglik + self.inside_bonus(sn)

            wSpotCell[:, n] = term_1 + term_2 + logeta_bar + loglik

        logger.info('end loop')
        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.parent_cell_prob = pSpotNeighb
        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()
        assert np.isfinite(wSpotCell).all(), "wSpotCell array contains non numeric values"
        logger.info(' Spot to cell loglikelihood: %f' % wSpotCell.sum())
        # self.spots.update_cell_prob(pSpotNeighb, self.cells)
        # logger.info('spot ---> cell probabilities updated')

    def spots_to_cell_par(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        # nN = self.config['nNeighbors'] + 1
        nN = self.nN
        nS = self.spots.data.gene_name.shape[0]

        misread_density_adj_factor = -1 * np.log(2 * np.pi * self.cells.mcr ** 2) / 2
        misread_density_adj = np.exp(misread_density_adj_factor) * self.config['MisreadDensity']
        # initialise array with the misread density
        wSpotCell = np.zeros([nS, nN]) + np.log(misread_density_adj)
        gn = self.spots.data.gene_name.values
        sc_mean = self.single_cell.log_mean_expression.loc[gn].values
        logeta_bar = self.genes.logeta_bar[self.spots.gene_id]

        # logger.info('start multiprocessing')
        # logger.info('start zip')
        args_in = list(zip(np.arange(nN - 1), [logeta_bar] * len(np.arange(nN - 1)), [sc_mean] * len(np.arange(nN - 1))))
        # logger.info('end zip')
        pool = ThreadPool()
        pout = pool.map(self.calc_loglik, args_in)
        pool.close()
        pool.join()
        for i, d in enumerate(pout):
            assert i == d[0]
            wSpotCell[:, i] = d[1]
        # logger.info('end multiprocessing')

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.parent_cell_prob = pSpotNeighb
        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()
        assert np.isfinite(wSpotCell).all(), "wSpotCell array contains non numeric values"
        logger.info(' Spot to cell loglikelihood: %f' % wSpotCell.sum())
        # self.spots.update_cell_prob(pSpotNeighb, self.cells)
        # logger.info('spot ---> cell probabilities updated')

    def inside_bonus(self, cell_id):
        mask = self.spots.data.label.values == cell_id
        bonus = self.config['InsideCellBonus']
        return np.where(mask, bonus, 0)

    def calc_loglik(self, args):
        n = args[0]
        logeta_bar = args[1]
        sc_mean = args[2]

        # get the spots' nth-closest cell
        sn = self.spots.parent_cell_id[:, n]

        # get the respective cell type probabilities
        cp = self.cells.classProb[sn]

        # multiply and sum over cells
        term_1 = np.einsum('ij, ij -> i', sc_mean, cp)

        log_gamma_bar = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

        term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)
        # logger.info('start mvn_loglik')
        loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells)

        # if 2D, add the bonus for the spots inside the boundaries.
        if not self.config['is_3D']:
            loglik = loglik + self.inside_bonus(sn)

        return n, term_1 + term_2 + logeta_bar + loglik

    # -------------------------------------------------------------------- #
    def eta_upd(self):
        """
        Calcs the expected eta
        Implements equation (5) of the Qian paper
        """
        grand_total = self.cells.background_counts.sum() + self.cells.total_counts.sum()
        assert round(grand_total) == self.spots.data.shape[0], \
            'The sum of the background spots and the total gene counts should be equal to the number of spots'

        classProb = self.cells.classProb
        mu = self.single_cell.mean_expression
        area_factor = self.cells.ini_cell_props['area_factor']
        gamma_bar = self.spots.gamma_bar

        zero_prob = classProb[:, -1]  # probability a cell being a zero expressing cell
        zero_class_counts = self.spots.zero_class_counts(self.spots.gene_id, zero_prob)

        # Calcs the sum in the Gamma distribution (equation 5). The zero class
        # is excluded from the sum, hence the arrays in the einsum below stop at :-1
        # Note. I also think we should exclude the "cell" that is meant to keep the
        # misreads, ie exclude the background
        class_total_counts = np.einsum('ck, gk, c, cgk -> g',
                                       classProb[:, :-1],
                                       mu.values[:, :-1],
                                       area_factor,
                                       gamma_bar[:, :, :-1])
        background_counts = self.cells.background_counts
        alpha = self.config['rGene'] + self.spots.counts_per_gene - background_counts - zero_class_counts
        beta = self.config['rGene'] / self.config['Inefficiency'] + class_total_counts

        # Finally, update gene efficiency
        self.genes.calc_eta(alpha, beta)

    # -------------------------------------------------------------------- #
    def gaussian_upd(self):
        self.centroid_upd()
        self.cov_upd()

    # -------------------------------------------------------------------- #
    def centroid_upd(self):
        spots = self.spots

        # get the total gene counts per cell
        N_c = self.cells.total_counts

        xyz_spots = spots.xyz_coords
        prob = spots.parent_cell_prob
        n = self.cells.config['nNeighbors'] + 1

        # mulitply the x coord of the spots by the cell prob
        a = np.tile(xyz_spots[:, 0], (n, 1)).T * prob

        # mulitply the y coord of the spots by the cell prob
        b = np.tile(xyz_spots[:, 1], (n, 1)).T * prob

        # mulitply the z coord of the spots by the cell prob
        c = np.tile(xyz_spots[:, 2], (n, 1)).T * prob

        # aggregated x and y coordinate
        idx = spots.parent_cell_id
        x_agg = npg.aggregate(idx.ravel(), a.ravel(), size=len(N_c))
        y_agg = npg.aggregate(idx.ravel(), b.ravel(), size=len(N_c))
        z_agg = npg.aggregate(idx.ravel(), c.ravel(), size=len(N_c))

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
        z_bar = np.nan * np.ones(N_c.shape)

        x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
        y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]
        z_bar[N_c > 0] = z_agg[N_c > 0] / N_c[N_c > 0]

        # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
        xyz_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T, z_bar.T)))

        # if you have a value for the estimated centroid use that, otherwise
        # use the initial (starting values) centroids
        ini_cent = self.cells.ini_centroids()
        xyz_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y'], ini_cent['z']])))

        # # sanity check. NaNs or Infs should appear together
        # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
        # use the fitted centroids where possible otherwise use the initial ones
        xyz_bar[np.isfinite(x_bar)] = xyz_bar_fitted[np.isfinite(x_bar)]
        self.cells.centroid = pd.DataFrame(xyz_bar, columns=['x', 'y', 'z'])
        # print(np.array(list(zip(x_bar.T, y_bar.T))))

    # -------------------------------------------------------------------- #
    def cov_upd(self):
        spots = self.spots

        # first get the scatter matrix
        S = self.cells.scatter_matrix(spots)  # sample sum of squares
        cov_0 = self.cells.ini_cov()
        nu_0 = self.cells.nu_0
        S_0 = cov_0 * nu_0  # prior sum of squarea
        N_c = self.cells.total_counts
        d = 2
        denom = np.ones([self.nC, 1, 1])
        denom[:, 0, 0] = N_c + nu_0
        # Note: need to add code to handle the case N_c + nu_0 <= d + 2
        cov = (S + S_0) / denom

        # delta = self.cells.ledoit_wolf(self.spots, cov)
        # stein = self.cells.stein(cov)
        # delta = self.cells.rblw(cov)

        delta = 0.5
        # logger.info('Mean shrinkage %4.2f' % delta.mean())
        # logger.info('cell 601 shrinkage %4.2f, %4.2f' % (shrinkage[601], sh[601]))
        #  logger.info('cell 601 gene counts %d' % self.cells.total_counts[601])
        # logger.info('cell 605 shrinkage %4.2f, %4.2f' % (shrinkage[605], sh[605]))
        #  logger.info('cell 605 gene counts %d' % self.cells.total_counts[605])
        # logger.info('cell 610 shrinkage %4.2f, %4.2f' % (shrinkage[610], sh[610]))
        #  logger.info('cell 610 gene counts %d' % self.cells.total_counts[610])

        # delta = delta.reshape(self.nC, 1, 1)
        # target = [np.trace(d)/2 * np.eye(2) for d in cov]  # shrinkage target
        cov = delta * cov_0 + (1 - delta) * cov
        # cov = delta * target + (1 - delta) * cov
        self.cells.cov = cov

    # -------------------------------------------------------------------- #
    def mu_upd(self):
        logger.info('Update single cell data')
        # make an array nS-by-nN and fill it with the spots id
        gene_ids = np.tile(self.spots.gene_id, (self.nN, 1)).T

        # flatten it
        gene_ids = gene_ids.ravel()

        # make corresponding arrays for cell_id and probs
        cell_ids = self.spots.parent_cell_id.ravel()
        probs = self.spots.parent_cell_prob.ravel()

        # make the array to be used as index in the group-by operation
        group_idx = np.vstack((cell_ids, gene_ids))

        # For each cell aggregate the number of spots from the same gene.
        # It will produce an array of size nC-by-nG where the entry at (c,g)
        # is the gene counts of gene g within cell c
        N_cg = npg.aggregate(group_idx, probs, size=(self.nC, self.nG))

        classProb = self.cells.classProb[1:, :-1].copy()
        geneCount = self.cells.geneCount[1:, :].copy()
        gamma_bar = self.spots.gamma_bar[1:, :, :-1].copy()
        area_factor = self.cells.ini_cell_props['area_factor'][1:]

        numer = np.einsum('ck, cg -> gk', classProb, geneCount)
        denom = np.einsum('ck, c, cgk, g -> gk', classProb, area_factor, gamma_bar, self.genes.eta_bar)

        # set the hyperparameter for the gamma prior
        # m = 1

        # ignore the last class, it is the zero class
        # numer = numer[:, :-1]
        # denom = denom[:, :-1]
        # mu_gk = (numer + m * self.single_cell.raw_data) / (denom + m)
        me, lme = self.single_cell._gene_expressions(numer, denom)
        self.single_cell._mean_expression = me
        self.single_cell._log_mean_expression = lme

        logger.info('Singe cell data updated')

    # -------------------------------------------------------------------- #
    def dalpha_upd(self):
        # logger.info('Update cell type (marginal) distribution')
        zeta = self.cells.classProb.sum(axis=0)  # this the class size
        alpha = self.cellTypes.ini_alpha()
        out = zeta + alpha

        # logger.info("***************************************************")
        # logger.info("**** Dirichlet alpha is zero if class size <=%d ****" % self.config['min_class_size'])
        # logger.info("***************************************************")
        mask = zeta <= self.config['min_class_size']
        out[mask] = 10e-6
        self.cellTypes.alpha = out

    # -------------------------------------------------------------------- #
    def redis_save(self):
        logger.info("redis start")
        # cell_prob_df = pd.DataFrame(data=self.spots.parent_cell_prob).set_index('spot_id')
        # self.redis_db.to_redis(cell_prob_df, "parent_cell_prob", iteration=self.iter_num,
        #                        has_converged=self.has_converged, unix_time=time.time())
        # self.redis_db.to_redis(self.spots.parent_cell_id, "parent_cell_id", iter_num=self.iter_num,
        #                        has_converged=self.has_converged, unix_time=time.time())
        # self.redis_db.to_redis(self.cells.geneCount, "geneCounts", iter_num=self.iter_num,
        #                        has_converged=self.has_converged, unix_time=time.time())
        # self.redis_db.to_redis(self.cells.classProb, "class_prob", iter_num=self.iter_num,
        #                        has_converged=self.has_converged, unix_time=time.time())

        eta_bar_df = pd.DataFrame({
            'gene_efficiency': self.genes.eta_bar,
            'gene': self.genes.gene_panel
        })
        self.redis_db.to_redis(eta_bar_df, "gene_efficiency", iteration=self.iter_num,
                               has_converged=self.has_converged, unix_time=time.time())

        pi_bar_df = pd.DataFrame({
            'weight': self.cellTypes.pi_bar,
            'class': self.cellTypes.names
        })
        self.redis_db.to_redis(pi_bar_df, "cell_type_prior", iter_num=self.iter_num,
                               has_converged=self.has_converged, unix_time=time.time())

        idx=[]
        for i, row in enumerate(self.cells.classProb):
            idx.append(np.argmax(row))
        prob = np.bincount(idx) / np.bincount(idx).sum()
        df = pd.DataFrame({
            'class_name': self.cellTypes.names,
            'prob': prob
        })
        self.redis_db.to_redis(df, "cell_type_posterior", iter_num=self.iter_num, has_converged=self.has_converged, unix_time=time.time())

        logger.info("redis end")
        # print(pd.DataFrame(from_redis(self.redis_db, "gene_efficiency")))

    def db_save(self):
        if self.conn is not None:
            try:
                self._db_save()
            except:
                self.conn.close()
                raise

    # -------------------------------------------------------------------- #
    def _db_save(self):
        db_opts = {'if_table_exists': 'replace'}  # choose between 'fail', 'replace', 'append'. Appending might make sense only if you want to see how estimates change from one loop to the next
        logger.info("starting db_save_geneCounts")
        self.db_save_geneCounts(self.iter_num, self.has_converged, db_opts)
        # self.db_save_cellByclass_prob(self.iter_num, self.has_converged, db_opts)
        logger.info("starting db_save_class_prob")
        self.db_save_class_prob(self.iter_num, self.has_converged, db_opts)
        logger.info("starting db_save_cell_type_posterior")
        self.db_save_cell_type_posterior(self.iter_num, self.has_converged, db_opts)
        # logger.info("starting db_save_parent_cell_prob")
        # self.db_save_parent_cell_prob(self.iter_num, self.has_converged, db_opts)
        # logger.info("starting db_save_parent_cell_id")
        # self.db_save_parent_cell_id(self.iter_num, self.has_converged, db_opts)
        logger.info("starting .genes.db_save")
        self.genes.db_save(self.conn, self.iter_num, self.has_converged, db_opts)
        logger.info("starting cellTypes.db_save")
        self.cellTypes.db_save(self.conn, self.iter_num, self.has_converged, db_opts)
        # self.spots.db_save(self.conn)
        # self.single_cell.db_save(self.conn, self.iter_num, self.has_converged, db_opts)

    # -------------------------------------------------------------------- #
    # def db_save_geneCounts(self, iter, has_converged, db_opts):
    #     df = pd.DataFrame(data=self.cells.geneCount,
    #                       index=np.arange(self.nC),
    #                       columns=self.genes.gene_panel)
    #     df.index.name = 'cell_label'
    #     df['iteration'] = iter
    #     df['has_converged'] = has_converged
    #     df['utc'] = datetime.datetime.utcnow()
    #     chunk_size = 999 // (len(df.columns) + 1)
    #     df.to_sql(name='geneCount', con=self.conn, if_exists=db_opts['if_table_exists'], chunksize=1000, method='multi')
    #     self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_label_iteration ON geneCount("cell_label", "iteration");')
    #
    # # -------------------------------------------------------------------- #
    # def db_save_class_prob(self, iter, has_converged, db_opts):
    #     df = pd.DataFrame(data=self.cells.classProb,
    #                       index=np.arange(self.nC),
    #                       columns=self.cellTypes.names)
    #     df.index.name = 'cell_label'
    #     df['iteration'] = iter
    #     df['has_converged'] = has_converged
    #     df['utc'] = datetime.datetime.utcnow()
    #     chunk_size = 999 // (len(df.columns) + 1)
    #     df.to_sql(name='class_prob', con=self.conn, if_exists=db_opts['if_table_exists'], chunksize=1000, method='multi')
    #     self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_label_iteration ON class_prob("cell_label", "iteration");')

    # -------------------------------------------------------------------- #
    def db_save_cell_type_posterior(self, iter, has_converged, db_opts):
        idx=[]
        for i, row in enumerate(self.cells.classProb):
            idx.append(np.argmax(row))
        prob = np.bincount(idx) / np.bincount(idx).sum()
        df = pd.DataFrame({
            'class_name': self.cellTypes.names,
            'prob': prob
        }).set_index('class_name')
        self.redis_db.to_redis(df, "cell_type_posterior", iteration=iter, has_converged=has_converged, unix_time=time.time())

    # -------------------------------------------------------------------- #
    # def db_save_parent_cell_prob(self, iter_num, has_converged, db_opts):
    #     df = pd.DataFrame(data=self.spots.parent_cell_prob,
    #                       index=np.arange(self.nS))
    #     df.index.name = 'spot_id'
    #     df['iteration'] = iter_num
    #     df['has_converged'] = has_converged
    #     df['utc'] = datetime.datetime.utcnow()
    #     chunk_size = 999 // (len(df.columns) + 1)
    #     df.to_sql(name='parent_cell_prob', con=self.conn, if_exists=db_opts['if_table_exists'], chunksize=1000, method='multi')
    #     self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_cell_ID_iteration ON parent_cell_prob("spot_id", "iteration");')
    #
    # # -------------------------------------------------------------------- #
    # def db_save_parent_cell_id(self, iter_num, has_converged, db_opts):
    #     df = pd.DataFrame(data=self.spots.parent_cell_id,
    #                       index=np.arange(self.nS))
    #     df.index.name = 'spot_id'
    #     df['iteration'] = iter_num
    #     df['has_converged'] = has_converged
    #     df['utc'] = datetime.datetime.utcnow()
    #     chunk_size = 999 // (len(df.columns) + 1)
    #     df.to_sql(name='parent_cell_id', con=self.conn, if_exists=db_opts['if_table_exists'], chunksize=1000, method='multi')
    #     self.conn.execute(
    #         'CREATE UNIQUE INDEX IF NOT EXISTS ix_cell_ID_iteration ON parent_cell_id("spot_id", "iteration");')

    def counts_within_radius(self, r):
        """
        calcs the gene counts within radius r of the centroid of each cell.
        Radius is in micron
        """

        # check that the background (label=0) is at the top row
        assert self.cells.ini_cell_props['cell_label'][0] == 0

        # radius r is assumed to be in micron. Convert it to pixels
        r_px = r * self.config['anisotropy']

        spots = pd.concat([self.spots.data, self.spots.data_excluded])
        gene_names, gene_id = np.unique(spots.gene_name.values, return_inverse=True)  # that needs to be encapsulated!!! Make a function in the class to set the ids
        spots = spots.assign(gene_id=gene_id)

        xyz_coords = spots[['x', 'y', 'z']].values
        point_tree = spatial.cKDTree(xyz_coords)
        nearby_spots = point_tree.query_ball_point(self.cells.centroid, r_px)

        out = np.zeros([self.cells.centroid.shape[0], len(gene_names)])

        for i, d in enumerate(nearby_spots):
            t = spots.gene_id[d]
            b = np.bincount(t)
            if  len(gene_names) - len(b) < 0:
                print('oops')
            out[i, :] = np.pad(b, (0, len(gene_names) - len(b)), 'constant')

        # for each cell get the most likely cell type
        cell_type = []
        for i, d in enumerate(self.cells.classProb):
            j = self.cells.class_names[np.argmax(d)]
            cell_type.append(j)

        temp = pd.DataFrame({
            'cell_label': self.cells.ini_cell_props['cell_label'],
            'cell_label_old': self.cells.ini_cell_props['cell_label_old'],
            'cell_type': cell_type
        })

        # assert np.all(temp.cell_label==out.index)
        df = pd.DataFrame(out, index=self.cells.centroid.index, columns=gene_names)
        df = pd.merge(temp, df, right_index=True, left_on='cell_label', how='left')

        # ignore the first row, it is the background
        return df.iloc[1:, :]

