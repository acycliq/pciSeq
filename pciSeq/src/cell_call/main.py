import os
import numpy as np
import pandas as pd
import numpy_groupies as npg
import sqlite3
import datetime
from pciSeq.src.cell_call.datatypes import Cells, Spots, Genes, SingleCell, CellType
from pciSeq.src.cell_call.summary import collect_data
import pciSeq.src.cell_call.utils as utils
from pciSeq.src.cell_call.log_config import logger


class VarBayes(object):
    def __init__(self, _cells_df, _spots_df, scRNAseq, config):
        self.config = config
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
        self.nN = self.config['nNeighbors'] + 1     # number of closest nearby cells, candidates for being parent
                                                    # cell of any given spot. The last cell will be used for the
                                                    # misread spots. (ie cell at position nN is the background)

        self.conn = self.db_connect(':memory:')  # or use 'pciSeq.db' to create a db on the filesystem
        self.iter_num = None
        self.has_converged = False

    def initialise(self):
        self.cellTypes.ini_prior()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.genes.eta = np.ones(self.nG) * self.config['Inefficiency']
        self.spots.parent_cell_id, _ = self.spots.cells_nearby(self.cells)
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id,
                                                              self.config)  # assign a spot to a cell if it is within its cell boundaries
        self.spots.gamma_bar = np.ones([self.nC, self.nG, self.nK]).astype(self.config['dtype'])

    # -------------------------------------------------------------------- #
    def run(self):
        p0 = None
        iss_df = None
        gene_df = None
        max_iter = self.config['max_iter']

        self.initialise()
        for i in range(max_iter):
            self.iter_num = i

            # 1. For each cell, calc the expected gene counts
            self.geneCount_upd()

            # 2. calc expected gamma
            self.gamma_upd()

            if self.config['is_3D']:
                self.gaussian_upd()

            # 3. assign cells to cell types
            self.cell_to_cellType()

            # update cell type weights
            if self.config['is_3D']:
                self.dalpha_upd()

            # 4. assign spots to cells
            self.spots_to_cell()

            # 5. update gene efficiency
            self.eta_upd()

            # 6. Update single cell data
            if self.single_cell.isMissing:
                self.mu_upd()

            self.has_converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info(' Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            if self.has_converged:
                self.db_save()
                iss_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell)
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
        beta = np.einsum('c, gk, g -> cgk', cells.cell_props['area_factor'], self.single_cell.mean_expression,
                         self.genes.eta).astype(dtype) + cfg['rSpot']
        # beta = np.einsum('c, gk -> cgk', cells.cell_props['area_factor'], self.single_cell.mean_expression).astype(dtype) + cfg['rSpot']
        rho = cfg['rSpot'] + cells.geneCount

        self.spots.gamma_bar = self.spots.gammaExpectation(rho, beta)
        self.spots.log_gamma_bar = self.spots.logGammaExpectation(rho, beta)

    # -------------------------------------------------------------------- #
    def cell_to_cellType(self):
        """
        return a an array of size numCells-by-numCellTypes where element in position [i,j]
        keeps the probability that cell i has cell type j
        Implements equation (2) of the Qian paper
        :param spots:
        :param config:
        :return:
        """

        # gene_gamma = self.genes.eta
        dtype = self.config['dtype']
        # ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.alpha, self.genes.eta, sc.mean_expression.data) + self.config['SpotReg']
        ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.cell_props['area_factor'], self.genes.eta,
                              self.single_cell.mean_expression).astype(dtype)
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
    def spots_to_cell_2D(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        nN = self.nN
        nS = self.spots.data.Gene.shape[0]

        # initialise array with the misread density
        aSpotCell = np.zeros([nS, nN])
        gn = self.spots.data.Gene.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values

        ## DN: 22-Jun-2022. I think this is missing from equation 4, Xiaoyan's paper
        ## I think we should multiply mu (single cell expression averages) by the gene efficiency
        unames, idx = np.unique(gn, return_inverse=True)
        assert np.all(unames == self.genes.gene_panel)

        eta = self.genes.eta[idx]
        expected_counts = np.einsum('sc, s -> sc', expected_counts,
                                    eta)  # multiply the single cell averages by the gene efficiency
        ## DN ends

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        my_D = np.zeros([nS, nN])
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            term_1 = np.einsum('ij, ij -> i', expected_counts, cp)

            # logger.info('genes.spotNo should be something like spots.geneNo instead!!')
            expectedLog = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, expectedLog)
            aSpotCell[:, n] = term_1 + term_2
            # my_D[:, n] = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells)
            my_covs = self.cells.cov[sn] * np.diag([1, 1, 1])

            # my_D[:, n] = self.spots.multiple_logpdfs(self.spots.xyz_coords[:, :2],  self.cells.centroid.values[sn][:, :2], my_covs[:,:2,:2])
            my_D[:, n] = self.spots.multiple_logpdfs(self.spots.xyz_coords, self.cells.centroid.values[sn], my_covs)
        wSpotCell = aSpotCell + self.spots.loglik(self.cells, self.config)

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.parent_cell_prob = pSpotNeighb
        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()
        # self.spots.update_cell_prob(pSpotNeighb, self.cells)
        # logger.info('spot ---> cell probabilities updated')

    # -------------------------------------------------------------------- #
    def spots_to_cell(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        # nN = self.config['nNeighbors'] + 1
        nN = self.nN
        nS = self.spots.data.Gene.shape[0]

        misread_density_adj_factor = -1 * np.log(2 * np.pi * self.cells.mcr ** 2) / 2
        misread_density_adj = np.exp(misread_density_adj_factor) * self.config['MisreadDensity']
        # initialise array with the misread density
        wSpotCell = np.zeros([nS, nN]) + np.log(misread_density_adj)
        gn = self.spots.data.Gene.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values

        ## DN: 22-Jun-2022. I think this is missing from equation 4, Xiaoyan's paper
        ## I think we should multiply mu (single cell expression averages) by the gene efficiency
        unames, idx = np.unique(gn, return_inverse=True)
        assert np.all(unames == self.genes.gene_panel)

        eta = self.genes.eta[idx]
        expected_counts = np.einsum('sc, s -> sc', expected_counts,
                                    eta)  # multiply the single cell averages by the gene efficiency
        ## DN ends

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        my_D = np.zeros([nS, nN])
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            term_1 = np.einsum('ij, ij -> i', expected_counts, cp)

            # logger.info('genes.spotNo should be something like spots.geneNo instead!!')
            expectedLog = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]
            # expectedLog = utils.bi2(self.elgamma.data, [nS, nK], sn[:, None], self.spots.data.gene_id.values[:, None])

            term_2 = np.einsum('ij, ij -> i', cp,
                               expectedLog)  # same as np.sum(cp * expectedLog, axis=1) but bit faster

            loglik = self.spots.mvn_loglik(self.spots.xyz_coords, sn, self.cells)
            my_D[:, n] = loglik
            wSpotCell[:, n] = term_1 + term_2 + loglik

        # update the prob a spot belongs to a neighboring cell
        pSpotNeighb = utils.softmax(wSpotCell, axis=1)
        self.spots.parent_cell_prob = pSpotNeighb
        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()
        assert np.isfinite(wSpotCell).all(), "wSpotCell array contains non numeric values"
        logger.info(' Spot to cell loglikelihood: %f' % wSpotCell.sum())
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

        classProb = self.cells.classProb
        mu = self.single_cell.mean_expression
        area_factor = self.cells.cell_props['area_factor']
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
        numer = self.config['rGene'] + self.spots.counts_per_gene - background_counts - zero_class_counts
        denom = self.config['rGene'] / self.config['Inefficiency'] + class_total_counts
        res = numer / denom

        # Finally, update gene efficiency
        self.genes.eta = res

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
        area_factor = self.cells.cell_props['area_factor'][1:]

        numer = np.einsum('ck, cg -> gk', classProb, geneCount)
        denom = np.einsum('ck, c, cgk, g -> gk', classProb, area_factor, gamma_bar, self.genes.eta)

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
    def db_save(self):
        try:
            self._db_save()
        except:
            self.conn.close()
            raise

    # -------------------------------------------------------------------- #
    def _db_save(self):
        db_opts = {'if_table_exists': 'replace'}  # choose between 'fail', 'replace', 'append'. Appending might make sense only if you want to see how estimates change from one loop to the next
        self.db_save_geneCounts(self.iter_num, self.has_converged, db_opts)
        self.db_save_class_prob(self.iter_num, self.has_converged, db_opts)
        self.db_save_parent_cell_prob(self.iter_num, self.has_converged, db_opts)
        self.db_save_parent_cell_id(self.iter_num, self.has_converged, db_opts)
        self.genes.db_save(self.conn, self.iter_num, self.has_converged, db_opts)
        self.cellTypes.db_save(self.conn, self.iter_num, self.has_converged, db_opts)
        self.spots.db_save(self.conn)
        self.single_cell.db_save(self.conn, self.iter_num, self.has_converged, db_opts)

    # -------------------------------------------------------------------- #
    def db_save_geneCounts(self, iter, has_converged, db_opts):
        df = pd.DataFrame(data=self.cells.geneCount,
                              index=np.arange(self.nC),
                              columns=self.genes.gene_panel)
        df.index.name = 'cell_label'
        df['iteration'] = iter
        df['has_converged'] = has_converged
        df['utc'] = datetime.datetime.utcnow()
        df.to_sql(name='geneCount', con=self.conn, if_exists=db_opts['if_table_exists'])
        self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_label_iteration ON geneCount("cell_label", "iteration");')

    # -------------------------------------------------------------------- #
    def db_save_class_prob(self, iter, has_converged, db_opts):
        df = pd.DataFrame(data=self.cells.classProb,
                          index=np.arange(self.nC),
                          columns=self.cellTypes.names)
        df.index.name = 'cell_label'
        df['iteration'] = iter
        df['has_converged'] = has_converged
        df['utc'] = datetime.datetime.utcnow()
        df.to_sql(name='classProb', con=self.conn, if_exists=db_opts['if_table_exists'])
        self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_label_iteration ON classProb("cell_label", "iteration");')

    # -------------------------------------------------------------------- #
    def db_save_parent_cell_prob(self, iter_num, has_converged, db_opts):
        df = pd.DataFrame(data=self.spots.parent_cell_prob,
                          index=np.arange(self.nS))
        df.index.name = 'spot_id'
        df['iteration'] = iter_num
        df['has_converged'] = has_converged
        df['utc'] = datetime.datetime.utcnow()
        df.to_sql(name='parent_cell_prob', con=self.conn, if_exists=db_opts['if_table_exists'])
        self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_cell_ID_iteration ON parent_cell_prob("spot_id", "iteration");')

    # -------------------------------------------------------------------- #
    def db_save_parent_cell_id(self, iter_num, has_converged, db_opts):
        df = pd.DataFrame(data=self.spots.parent_cell_id,
                          index=np.arange(self.nS))
        df.index.name = 'spot_id'
        df['iteration'] = iter_num
        df['has_converged'] = has_converged
        df['utc'] = datetime.datetime.utcnow()
        df.to_sql(name='parent_cell_id', con=self.conn, if_exists=db_opts['if_table_exists'])
        self.conn.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_cell_ID_iteration ON parent_cell_id("spot_id", "iteration");')

    # -------------------------------------------------------------------- #
    def db_connect(self, dbpath):
        if os.path.isfile(dbpath):
            os.remove(dbpath)
        return sqlite3.connect(dbpath)
