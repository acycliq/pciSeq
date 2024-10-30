import time
import numpy as np
import pandas as pd
import numpy_groupies as npg
import scipy.spatial as spatial
import json
import os
import time
from scipy.special import softmax
import pciSeq.src.core.utils as utils
from pciSeq.src.core.summary import collect_data
from pciSeq.src.diagnostics.utils import redis_db
from dask.delayed import delayed
import shutil
from pciSeq.src.core.datatypes import Cells, Spots, Genes, SingleCell, CellType
import logging

main_logger = logging.getLogger(__name__)


class VarBayes:
    def __init__(self, _cells_df, _spots_df, scRNAseq, config):
        self.config = config
        self.redis_db = redis_db(flush=True) if config['is_redis_running'] else None
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
        self.nN = self.config['nNeighbors'] + 1  # number of closest nearby cells, candidates for being parent
        # cell of any given spot. The last cell will be used for the
        # misread spots. (ie cell at position nN is the background)
        self.iter_num = None
        self.iter_delta = []
        self.has_converged = False
        self._scaled_exp = None
        self.needs_recalc = True

    @property
    def scaled_exp(self):
        return self._scaled_exp

    def __getstate__(self):
        # set here attributes to be excluded from serialisation (pickling)
        # It makes the pickle filesize smaller but maybe this will have to
        # change in the future.
        # Removing redis so I can pickle and _scaled_exp because it is delayed
        # FYI: https://realpython.com/python-pickle-module/
        attributes = self.__dict__.copy()
        del attributes['redis_db']
        # del attributes['_scaled_exp']
        return attributes

    def initialise(self):
        self.cellTypes.ini_prior()
        self.cells.classProb = np.tile(self.cellTypes.prior, (self.nC, 1))
        self.genes.init_eta(1, 1 / self.config['Inefficiency'])
        self.spots.parent_cell_id = self.spots.cells_nearby(self.cells)
        self.spots.parent_cell_prob = self.spots.ini_cellProb(self.spots.parent_cell_id, self.config)

    # -------------------------------------------------------------------- #
    def run(self):
        self.initialise()
        cell_df, gene_df = self.main_loop()
        return cell_df, gene_df

    def main_loop(self):
        p0 = None
        cell_df = None
        gene_df = None
        max_iter = self.config['max_iter']

        for i in range(max_iter):
            self.iter_num = i

            # 1. For each cell, calc the expected gene counts
            self.geneCount_upd()

            # 2. calc expected gamma
            self.gamma_upd()

            # 3 update correlation matrix and variance of the gaussian distribution
            if self.single_cell.isMissing or (self.config['InsideCellBonus'] is False):
                self.gaussian_upd()

            # 4. assign cells to cell types
            self.cell_to_cellType()

            # 5. assign spots to cells
            self.spots_to_cell()

            # 6. update gene efficiency
            # if 'overrides' not in self.config:
            if self.needs_recalc:
                self.eta_upd()

            # 7. update the dirichlet distribution
            if self.single_cell.isMissing or (self.config['cell_type_prior'] == 'weighted'):
                self.dalpha_upd()

            # 8. Update single cell data
            if self.single_cell.isMissing:
                self.mu_upd()

            self.has_converged, delta = utils.hasConverged(self.spots, p0, self.config['CellCallTolerance'])
            main_logger.info('Iteration %d, mean  prob change %f' % (i, delta))

            # keep track of the deltas
            self.iter_delta.append(delta)

            # replace p0 with the latest probabilities
            p0 = self.spots.parent_cell_prob

            if self.config['is_redis_running']:
                self.redis_upd()

            if self.has_converged:
                self.cell_analysis(2259)
                cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell)
                break

            if i == max_iter - 1:
                main_logger.info('Loop exhausted. Exiting with convergence status: %s' % self.has_converged)
        return cell_df, gene_df

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
        out = np.zeros([self.nC, self.nG], dtype=np.float32)
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

        self._scaled_exp = delayed(utils.scaled_exp(cells.ini_cell_props['area_factor'],
                                                    self.single_cell.mean_expression.values,
                                                    self.genes.eta_bar))

        beta = self.scaled_exp.compute() + cfg['rSpot']
        rho = cfg['rSpot'] + cells.geneCount

        self.spots._gamma_bar = delayed(self.spots.gammaExpectation(rho, beta))
        self.spots._log_gamma_bar = delayed(self.spots.logGammaExpectation(rho, beta))

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

        ScaledExp = self.scaled_exp.compute()
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)
        wCellClass = np.sum(contr, axis=1) + self.cellTypes.log_prior
        pCellClass = softmax(wCellClass, axis=1)
        del contr

        self.cells.classProb = pCellClass

    # -------------------------------------------------------------------- #
    def spots_to_cell(self):
        """
        spot to cell assignment.
        Implements equation (4) of the Qian paper
        """
        nN = self.nN
        nS = self.spots.data.gene_name.shape[0]

        wSpotCell = np.zeros([nS, nN], dtype=np.float32)
        gn = self.spots.data.gene_name.values
        expected_counts = self.single_cell.log_mean_expression.loc[gn].values
        logeta_bar = self.genes.logeta_bar[self.spots.gene_id]

        # loglik = self.spots.loglik(self.cells, self.config)

        # pre-populate last column
        wSpotCell[:, -1] = np.log(self.config['MisreadDensity'])

        # loop over the first nN-1 closest cells. The nN-th column is reserved for the misreads
        for n in range(nN - 1):
            # get the spots' nth-closest cell
            sn = self.spots.parent_cell_id[:, n]

            # get the respective cell type probabilities
            cp = self.cells.classProb[sn]

            # multiply and sum over cells
            term_1 = np.einsum('ij, ij -> i', expected_counts, cp)

            log_gamma_bar = self.spots.log_gamma_bar.compute()
            log_gamma_bar = log_gamma_bar[self.spots.parent_cell_id[:, n], self.spots.gene_id]

            term_2 = np.einsum('ij, ij -> i', cp, log_gamma_bar)

            # wSpotCell[:, n] = term_1 + term_2 + logeta_bar + loglik[:, n]
            mvn_loglik = self.spots.mvn_loglik(self.spots.xy_coords, sn, self.cells)
            wSpotCell[:, n] = term_1 + term_2 + logeta_bar + mvn_loglik

        # apply inside cell bonus
        # NOTE. This is not applied 100% correctly. For example the fourth spot in the demo data. Id2, (x, y) = (0, 4484)
        # The spots is within the boundaries of cell with label 5 but this is not its closest cell. The closest cell is
        # cell label = 4. Cell label=5 is the second closest cell and label=4 the first closest. Therefore, the bonus should be
        # applied when we handle the column for the seconds closest near-by cell. The implementation below implies that if
        # a spot is inside the cell boundaries then that cell is the closest one.
        mask = np.greater(self.spots.data.label.values, 0, where=~np.isnan(self.spots.data.label.values))
        wSpotCell[mask, 0] = wSpotCell[mask, 0] + self.config['InsideCellBonus']

        # update the prob a spot belongs to a neighboring cell
        self.spots.parent_cell_prob = softmax(wSpotCell, axis=1)

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
        gamma_bar = self.spots.gamma_bar.compute()

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

        # Finally, update gene_gamma
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

        # shrinkage
        delta = 0.5
        cov = delta * cov_0 + (1 - delta) * cov
        self.cells.cov = cov.astype(np.float32)

    # -------------------------------------------------------------------- #
    def mu_upd(self):
        classProb = self.cells.classProb[1:, :-1].copy()
        geneCount = self.cells.geneCount[1:, :].copy()
        gamma_bar = self.spots.gamma_bar.compute()[1:, :, :-1]
        area_factor = self.cells.ini_cell_props['area_factor'][1:]

        numer = np.einsum('ck, cg -> gk', classProb, geneCount)
        denom = np.einsum('ck, c, cgk, g -> gk', classProb, area_factor, gamma_bar, self.genes.eta_bar)

        me, lme = self.single_cell._gene_expressions(numer, denom)
        self.single_cell._mean_expression = me
        self.single_cell._log_mean_expression = lme

    # -------------------------------------------------------------------- #
    def dalpha_upd(self):
        # logger.info('Update cell type (marginal) distribution')
        zeta = self.cells.classProb.sum(axis=0)  # this the class size
        alpha = self.cellTypes.ini_alpha()
        out = zeta + alpha

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
        gene_names, gene_id = np.unique(self.spots.data.gene_name.values,
                                        return_inverse=True)  # that needs to be encapsulated!!! Make a function in the class to set the ids
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
        self.redis_db.publish(eta_bar_df, "gene_efficiency", iteration=self.iter_num,
                              has_converged=self.has_converged, unix_time=time.time())

        idx = []
        size = self.cellTypes.names.shape[0]
        for i, row in enumerate(self.cells.classProb[1:,
                                :]):  # ignore the top row, it corresponds to the background, it is not an actual cell
            idx.append(np.argmax(row))
        counts = np.bincount(idx, minlength=size)
        # prob = np.bincount(idx) / np.bincount(idx).sum()
        df = pd.DataFrame({
            'class_name': self.cellTypes.names,
            'counts': counts
        })
        self.redis_db.publish(df, "cell_type_posterior", iteration=self.iter_num, has_converged=self.has_converged,
                              unix_time=time.time())

    def calculate_spot_contributions(self, cell_num, datapoints):
        self.config['launch_diagnostics'] = False
        self.config['launch_viewer'] = False
        self.config['save_data'] = False
        self.config['is_redis_running'] = False
        self.needs_recalc = False

        df_before = self.get_celltypes_for_cell(cell_num)

        mask = np.any(self.spots.cells_nearby(self.cells) == cell_num, axis=1)
        df = self.spots.data[mask]

        self.spots.Dist = self.spots.Dist[mask]
        self.spots.data = self.spots.data[mask]
        self.spots.gene_id = self.spots.gene_id[mask]
        self.spots.parent_cell_id = self.spots.parent_cell_id[mask]
        self.spots.parent_cell_prob = self.spots.parent_cell_prob[mask]
        # self.spots.xy_coords = self.spots.xy_coords[mask]
        _, _, counts_per_gene_masked = np.unique(self.spots.data.gene_name.values, return_inverse=True,
                                                 return_counts=True)
        self.nS = self.spots.data.shape[0]
        self.spots.nS = self.spots.data.shape[0]

        datapoints = utils.get_closest(df, datapoints)
        if datapoints is None:
            df_after = df_before.copy()
        else:
            # set the overrides
            self.set_overrides(datapoints)
            self.spots.parent_cell_prob = self.spots.parent_cell_prob

            # redo the estimation
            _ = self.main_loop()

            df_after = self.get_celltypes_for_cell(cell_num)

        return df_before, df_after

    def set_overrides(self, data):
        self.config['overrides'] = data.to_dict(orient='records')
        main_logger.info('overrides saved!')

    def get_celltypes_for_cell(self, cell_num):
        cell_df, gene_df = collect_data(self.cells, self.spots, self.genes, self.single_cell)
        df = cell_df[cell_df.Cell_Num == cell_num]
        class_name = df.ClassName.tolist()[0]
        prob = df.Prob.tolist()[0]
        df = pd.DataFrame({
            'class_name': class_name,
            'prob': prob
        })
        return df

    def gene_loglik_contributions(self, cell_num, user_class=None):
        """
        Calculate and return gene log-likelihood contributions for a specified cell.

        This function analyzes how each gene contributes to the cell type classification
        by calculating log-likelihood values for each gene under different cell type hypotheses.

        Args:
            cell_num (int): The cell number to analyze. Must be between 0 and nC-1.
            user_class (str, optional): The cell class to compare against the assigned class.
                If None, uses the assigned class.

        Returns:
            dict: A dictionary containing:
                - assigned_class (str): The automatically assigned cell class
                - user_class (str): The user-specified class for comparison
                - assigned_contr (list): Log-likelihood contributions for assigned class
                - cell_num (int): The analyzed cell number
                - gene_names (list): List of gene names
                - class_names (list): List of available cell type classes
                - class_probs (dict): Probability distribution over cell types
                - contr (dict): Log-likelihood contributions for all classes

        Raises:
            ValueError: If cell_num is invalid or user_class is not recognized
            """

        if cell_num < 0 or cell_num >= self.nC:
            raise ValueError(f"Invalid cell number. Must be between 0 and {self.nC - 1}")

        assigned_class_idx = np.argmax(self.cells.classProb[cell_num])
        assigned_class = self.cellTypes.names[assigned_class_idx]

        if user_class is None:
            user_class = assigned_class

        try:
            user_class_idx = np.where(self.cellTypes.names == user_class)[0][0]
        except IndexError:
            raise ValueError(
                f"Invalid user class: {user_class}. Available classes are: {', '.join(self.cellTypes.names)}")

        # it will be nice to avoid duplicating the code. This has already been calculated.
        ScaledExp = self.scaled_exp.compute()
        pNegBin = ScaledExp / (self.config['rSpot'] + ScaledExp)
        cgc = self.cells.geneCount
        contr = utils.negBinLoglik(cgc, self.config['rSpot'], pNegBin)

        # Calculate contributions for all classes
        all_class_contrs = contr[cell_num, :, :]

        # Prepare the user_data dictionary with contributions for all classes
        user_data = {
            class_name: all_class_contrs[:, class_idx].tolist()
            for class_idx, class_name in enumerate(self.cellTypes.names)
        }

        class_probs = dict(zip(self.cellTypes.names.tolist(), self.cells.classProb[cell_num].tolist()))

        out = {
            'assigned_class': assigned_class,
            'user_class': user_class,
            'assigned_contr': all_class_contrs[:, assigned_class_idx].tolist(),
            'cell_num': cell_num,
            'gene_names': self.genes.gene_panel.tolist(),
            'class_names': self.cellTypes.names.tolist(),
            'class_probs': class_probs,
            'contr': user_data
        }

        # Call the plotting function
        # utils.gene_loglik_contributions_scatter(out)

        return out

    def spot_dist_and_prob(self, cell_num):
        """
        Calculate the relationship between spot-to-cell distances and their assignment
        probabilities for a given cell.

        This function analyzes spatial relationships between spots and a target cell by:
        1. Finding spots near the target cell
        2. Calculating distances from spots to cell centroid
        3. Computing assignment probabilities

        Args:
            cell_num (int): The cell number to analyze

        Returns:
            dict: A dictionary containing plot data:
                - x (list): Distances from spots to cell centroid
                - y (list): Assignment probabilities
                - labels (list): Gene names for each spot
                - cell_num (int): The analyzed cell number
                - title (str): Plot title
                - xlabel (str): X-axis label
                - ylabel (str): Y-axis label

        Note:
            The returned data is structured for visualization in the cell analysis dashboard.
        """

        # Get cell centroid
        centroid_yx = self.cells.yx_coords[cell_num]

        # Find spots near target cell
        is_spot_near_target_cell = self.spots.parent_cell_id == cell_num
        mask = np.any(is_spot_near_target_cell, axis=1)

        # Get probabilities for spots assigned to this cell
        prob = self.spots.parent_cell_prob[is_spot_near_target_cell]

        # Select relevant spots
        spots = self.spots.data[mask]

        # Find most likely parent cell for each spot
        max_idx = np.argmax(self.spots.parent_cell_prob[mask], axis=1)
        assigned_cell = np.choose(max_idx, self.spots.parent_cell_id[mask].T)

        # Calculate distances and create DataFrame
        spots = spots.assign(
            cell_x=centroid_yx[1],
            cell_y=centroid_yx[0],
            prob=prob,
            assigned_cell=assigned_cell,
            dist=np.sqrt((spots.x - centroid_yx[1]) ** 2 + (spots.y - centroid_yx[0]) ** 2),
        )

        # Select and order columns
        spots = spots[
            ['x', 'y', 'gene_name', 'cell_x', 'cell_y',
             'prob', 'assigned_cell', 'dist']
        ].reset_index()

        # Prepare plot data
        # Prepare plot data with cell-specific axis labels
        data = {
            'x': spots.dist.tolist(),
            'y': spots.prob.tolist(),
            'labels': spots.gene_name.tolist(),
            'cell_num': cell_num,
            'title': f'Cell {cell_num} - Distance vs Assignment Probability',
            'xlabel': f'Distance from cell {cell_num} centroid',
            'ylabel': f'Assignment probability to cell {cell_num}'
        }
        return data

    def cell_analysis(self, cell_num, output_dir=None):
        """
        Generates data and launches the cell analysis dashboard for a specific cell.

        Args:
            cell_num: The cell number to analyze
            output_dir: Optional directory to save JSON files. Defaults to cell_analysis/
        """
        # Get default output directory if none specified
        if output_dir is None:
            output_dir = utils.get_out_dir(self.config['output_path'])
            output_dir = os.path.join(output_dir, 'debug', 'cell_analysis')

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Get the data2
        assigned_class_idx = np.argmax(self.cells.classProb[cell_num])
        user_class = self.cellTypes.names[assigned_class_idx]

        # Generate gene contribution data
        spot_dist = self.spot_dist_and_prob(cell_num)

        loglik_data = self.gene_loglik_contributions(cell_num, user_class)
        with open(os.path.join(output_dir, 'gene_loglik_contr.json'), 'w') as fp:
            json.dump(loglik_data, fp)
            main_logger.info(f'saved at {os.path.join(output_dir, "gene_loglik_contr.json")}')

        # Save the data files
        with open(os.path.join(output_dir, "spot_dist.json"), "w") as f:
            json.dump(spot_dist, f)
            main_logger.info(f'saved at {os.path.join(output_dir, "spot_dist.json")}')

        pciSeq_dir = utils.get_pciSeq_install_dir()
        src = os.path.join(pciSeq_dir, 'static', 'cell_analysis')

        shutil.copytree(src, output_dir, dirs_exist_ok=True)
        # viewer_utils_logger.info('viewer code (%s) copied from %s to %s' % (dim, src, dst))

        # Launch the dashboard
        import webbrowser
        import http.server
        import socketserver
        import threading
        import random

        # Start HTTP server
        os.chdir(output_dir)
        PORT = 8000 + random.randint(0, 999)

        Handler = http.server.SimpleHTTPRequestHandler

        def start_server():
            with socketserver.TCPServer(("", PORT), Handler) as httpd:
                print(f"Serving cell analysis dashboard at http://localhost:{PORT}")
                httpd.serve_forever()

        # Start server in a separate thread
        server_thread = threading.Thread(target=start_server, daemon=True)
        server_thread.start()

        # Open the dashboard in the default browser.
        # Add the timestamp as version number to prevent loading from the cache
        webbrowser.open(f'http://localhost:{PORT}/dashboard/cell_index.html?v={time.time()}')

        # try:
        #     input("Press Enter to stop the server and close the dashboard...")
        # except KeyboardInterrupt:
        #     print("\nShutting down the server...")




