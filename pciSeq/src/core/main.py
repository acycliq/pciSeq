import numpy as np
import pandas as pd
import scipy.spatial as spatial
import time
import numpy_groupies as npg
from pciSeq.src.core.datatypes import Cells, Spots, Genes, SingleCell, CellType
from pciSeq.src.core.summary import collect_data
import pciSeq.src.core.utils as utils
from pciSeq.src.diagnostics.utils import redis_db
from pciSeq.src.core.log_config import logger


class VarBayes:
    def __init__(self, _cells_df, _spots_df, scRNAseq, config):
        self.config = config
        self.redis_db = redis_db(flush=True) if config['is_redis_running'] else None
        self.cells = Cells(_cells_df, config)
        self.spots = Spots(_spots_df, config)
        self.genes = Genes(self.spots)
        self.single_cell = SingleCell(scRNAseq, self.genes.gene_panel, self.config)
        self.cellTypes = CellType(self.single_cell)
        self.cells.class_names = self.single_cell.classes
        self.nC = self.cells.nC                         # number of cells
        self.nG = self.genes.nG                         # number of genes
        self.nK = self.cellTypes.nK                     # number of classes
        self.nS = self.spots.nS                         # number of spots
        self.nN = self.config['nNeighbors'] + 1         # number of closest nearby cells, candidates for being parent
                                                        # cell of any given spot. The last cell will be used for the
                                                        # misread spots. (ie cell at position nN is the background)
        self.iter_num = None
        self.has_converged = False

    def __getstate__(self):
        # set here attributes to be excluded from serialisation (pickling)
        # It makes the pickle filesize smaller but maybe this will have to
        # change in the future.
        # These two attributes take up a lot of space on the disk:
        # _gamma_bar and _log_gamma_bar
        # FYI: https://realpython.com/python-pickle-module/
        attributes = self.__dict__.copy()
        del attributes['redis_db']
        return attributes

    def initialise(self):
        self.cellTypes.ini_prior()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.genes.init_eta(1, 1 / self.config['Inefficiency'])
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id, self.config)
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

            # 3. assign cells to cell types
            self.cell_to_cellType()

            # 4. assign spots to cells
            self.spots_to_cell()

            # 5. update gene efficiency
            self.eta_upd()

            # 6. Update single cell data
            if self.single_cell.isMissing:
                # 6.1 update correlation matrix and variance of the gaussian distribution
                self.gaussian_upd()
                # 6.2 update the dirichlet distribution
                self.dalpha_upd()
                # 6.3 update single cell data
                self.mu_upd()

            self.has_converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            logger.info(' Iteration %d, mean prob change %f' % (i, delta))

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            if self.config['is_redis_running']:
                self.redis_upd()

            if self.has_converged:
                # self.counts_within_radius(20)
                iss_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell)
                break

            if i == max_iter-1:
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
        beta = np.einsum('c, gk, g -> cgk', cells.ini_cell_props['area_factor'], self.single_cell.mean_expression, self.genes.eta_bar).astype(dtype) + cfg['rSpot']
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
        ScaledExp = np.einsum('c, g, gk -> cgk', self.cells.ini_cell_props['area_factor'], self.genes.eta_bar, self.single_cell.mean_expression).astype(dtype)
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cellTypes.log_prior
        pCellClass = utils.softmax(wCellClass, axis=1)

        ## self.cells.classProb = pCellClass
        # logger.info('Cell 0 is classified as %s with prob %4.8f' % (
        #     self.cells.class_names[np.argmax(wCellClass[0, :])], pCellClass[0, np.argmax(wCellClass[0, :])]))
        # logger.info('cell ---> cell class probabilities updated')
        self.cells.classProb = pCellClass
        # return pCellClass

    # -------------------------------------------------------------------- #
    def spots_to_cell(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        nN = self.nN
        nS = self.spots.data.gene_name.shape[0]

        aSpotCell = np.zeros([nS, nN])
        gn = self.spots.data.gene_name.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values
        logeta_bar = self.genes.logeta_bar[self.spots.gene_id]

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            term_1 = np.einsum('ij, ij -> i', expected_counts, cp)

            # logger.info('genes.spotNo should be something like spots.geneNo instead!!')
            log_gamma_bar = self.spots.log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)
            aSpotCell[:, n] = term_1 + term_2 + logeta_bar
        wSpotCell = aSpotCell + self.spots.loglik(self.cells, self.config)

        # update the prob a spot belongs to a neighboring cell
        self.spots.parent_cell_prob = utils.softmax(wSpotCell, axis=1)

        # Since the spot-to-cell assignments changed you need to update the gene counts now
        self.geneCount_upd()

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
        beta = self.config['rGene']/self.config['Inefficiency'] + class_total_counts
        # res = num / denom

        # Finally, update gene_gamma
        self.genes.calc_eta(alpha, beta)
        # self.genes.eta_bar = res.values

    # -------------------------------------------------------------------- #
    def gaussian_upd(self):
        self.centroid_upd()
        self.cov_upd()

    # -------------------------------------------------------------------- #
    def centroid_upd(self):
        spots = self.spots

        # get the total gene counts per cell
        N_c = self.cells.total_counts

        xy_spots = spots.xy_coords
        prob = spots.parent_cell_prob
        n = self.cells.config['nNeighbors'] + 1

        # multiply the x coord of the spots by the cell prob
        a = np.tile(xy_spots[:, 0], (n, 1)).T * prob

        # multiply the y coord of the spots by the cell prob
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
        S_0 = cov_0 * nu_0  # prior sum of squares
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
        # logger.info('Update single cell data')
        # # make an array nS-by-nN and fill it with the spots id
        # gene_ids = np.tile(self.spots.gene_id, (self.nN, 1)).T
        #
        # # flatten it
        # gene_ids = gene_ids.ravel()
        #
        # # make corresponding arrays for cell_id and probs
        # cell_ids = self.spots.parent_cell_id.ravel()
        # probs = self.spots.parent_cell_prob.ravel()
        #
        # # make the array to be used as index in the group-by operation
        # group_idx = np.vstack((cell_ids, gene_ids))
        #
        # # For each cell aggregate the number of spots from the same gene.
        # # It will produce an array of size nC-by-nG where the entry at (c,g)
        # # is the gene counts of gene g within cell c
        # N_cg = npg.aggregate(group_idx, probs, size=(self.nC, self.nG))

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

        # logger.info('Single cell data updated')

    # -------------------------------------------------------------------- #
    def dalpha_upd(self):
        # logger.info('Update cell type (marginal) distribution')
        zeta = self.cells.classProb.sum(axis=0)  # this the class size
        alpha = self.cellTypes.ini_alpha()
        out = zeta + alpha

        # logger.info("***************************************************")
        # logger.info("**** Dirichlet alpha is zero if class size <=%d ****" % self.config['min_class_size'])
        # logger.info("***************************************************")

        # 07-May-2023: Hiding 'min_class_size' from the config file. Should bring it back at a later version
        # mask = zeta <= self.config['min_class_size']
        min_class_size = 5
        mask = zeta <= min_class_size

        # make sure Zero class is the last one
        assert self.cellTypes.names[-1] == "Zero"
        assert len(self.cellTypes.names) == len(mask)

        # make sure the last value ie always False, overriding if necessary the
        # check a few lines above when the mask variable was set.
        # In this manner we will prevent the Zero class from being removed.
        mask[-1] = False

        # If a class size is smaller than 'min_class_size' then it will be assigned a weight of almost zero
        out[mask] = 10e-6
        self.cellTypes.alpha = out

    def counts_within_radius(self, r):
        """
        calcs the gene counts within radius r of the centroid of each cell.
        Units in radius are the same as in your coordinates x and y of your spots.
        """

        # check that the background (label=0) is at the top row
        assert self.cells.ini_cell_props['cell_label'][0] == 0

        # spots = pd.concat([self.spots.data, self.spots.data_excluded])
        gene_names, gene_id = np.unique(self.spots.data.gene_name.values, return_inverse=True)  # that needs to be encapsulated!!! Make a function in the class to set the ids
        spots = self.spots.data.assign(gene_id=gene_id)

        xy_coords = spots[['x', 'y']].values
        point_tree = spatial.cKDTree(xy_coords)
        nearby_spots = point_tree.query_ball_point(self.cells.centroid, r)

        out = np.zeros([self.cells.centroid.shape[0], len(gene_names)])

        for i, d in enumerate(nearby_spots):
            t = spots.gene_id[d]
            b = np.bincount(t)
            if len(gene_names) - len(b) < 0:
                print('oops')
            out[i, :] = np.pad(b, (0, len(gene_names) - len(b)), 'constant')

        # for each cell get the most likely cell type
        cell_type = []
        for i, d in enumerate(self.cells.classProb):
            j = self.cells.class_names[np.argmax(d)]
            cell_type.append(j)

        temp = pd.DataFrame({
            'cell_label': self.cells.ini_cell_props['cell_label'],
            'cell_type': cell_type
        })

        # if labels have been reassigned, there should be a column 'cell_label_old'.
        # Relabelling will happen if for example the original labelling that was
        # assigned from segmentation is not a continuous sequence of integers
        # Append the original labels to the temp folder so there is a mapping between
        # old and new labels
        if 'cell_label_old' in self.cells.ini_cell_props.keys():
            temp['cell_label_old'] = self.cells.ini_cell_props['cell_label_old']
            temp = temp[['cell_label', 'cell_label_old', 'cell_type']]

        # assert np.all(temp.cell_label==out.index)
        df = pd.DataFrame(out, index=self.cells.centroid.index, columns=gene_names)
        df = pd.merge(temp, df, right_index=True, left_on='cell_label', how='left')

        # ignore the first row, it is the background
        return df.iloc[1:, :]

    # -------------------------------------------------------------------- #
    def redis_upd(self):
        eta_bar_df = pd.DataFrame({
            'gene_efficiency': self.genes.eta_bar,
            'gene': self.genes.gene_panel
        })
        self.redis_db.to_redis(eta_bar_df, "gene_efficiency", iteration=self.iter_num,
                               has_converged=self.has_converged, unix_time=time.time())
        self.redis_db.publish(eta_bar_df, "gene_efficiency", iteration=self.iter_num,
                               has_converged=self.has_converged, unix_time=time.time())

        # pi_bar_df = pd.DataFrame({
        #     'weight': np.bincount(self.cells.classProb.argmax(axis=1))/np.bincount(self.cells.classProb.argmax(axis=1)).sum(),
        #     'class': self.cellTypes.names
        # })
        # self.redis_db.to_redis(pi_bar_df, "cell_type_prior", iter_num=self.iter_num,
        #                        has_converged=self.has_converged, unix_time=time.time())

        idx=[]
        size = self.cellTypes.names.shape[0]
        for i, row in enumerate(self.cells.classProb[1:, :]):  # ignore the top row, it corresponds to the background, it is not an actual cell
            idx.append(np.argmax(row))
        counts = np.bincount(idx, minlength=size)
        # prob = np.bincount(idx) / np.bincount(idx).sum()
        df = pd.DataFrame({
            'class_name': self.cellTypes.names,
            'counts': counts
        })
        self.redis_db.to_redis(df, "cell_type_posterior", iter_num=self.iter_num, has_converged=self.has_converged, unix_time=time.time())
        self.redis_db.publish(df, "cell_type_posterior", iteration=self.iter_num, has_converged=self.has_converged, unix_time=time.time())







