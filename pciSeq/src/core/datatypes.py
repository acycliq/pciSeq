import sys
import numpy as np
import pandas as pd
import scipy
import numexpr as ne
import numpy_groupies as npg
from natsort import natsort_keygen
from sklearn.neighbors import NearestNeighbors
from pciSeq.src.core.log_config import logger


class Cells(object):
    # Get rid of the properties where not necessary!!
    def __init__(self, _cells_df, config):
        self.config = config
        self.ini_cell_props, self.mcr = self.read_image_objects(_cells_df, config)
        self.nC = len(self.ini_cell_props['cell_label'])
        self.classProb = None
        self.class_names = None
        self._cov = self.ini_cov()
        self.nu_0 = config['mean_gene_counts_per_cell']
        self._centroid = self.ini_centroids()
        self._gene_counts = None
        self._background_counts = None

    # -------- PROPERTIES -------- #
    @property
    def yx_coords(self):
        return self.centroid[['y', 'x']].values

    @property
    def geneCount(self):
        return self._gene_counts

    @geneCount.setter
    def geneCount(self, val):
        self._gene_counts = val

    @property
    def background_counts(self):
        return self._background_counts

    @background_counts.setter
    def background_counts(self, val):
        # assert val[1:, :].sum() == 0, 'Input array must be zero everywhere apart from the top row'
        # self._background_counts = val[0, :]
        self._background_counts = val

    @property
    def total_counts(self):
        # tc = self.geneCount.sum(axis=1)
        return self.geneCount.sum(axis=1)

    @property
    def centroid(self):
        # lst = list(zip(*[self._centroid['x'], self._centroid['y']]))
        return self._centroid.copy()

    @centroid.setter
    def centroid(self, df):
        assert isinstance(df, pd.DataFrame), 'Input should be a dataframe'
        assert set(df.columns.values) == {'x', 'y'}, 'Dataframe columns should be ''x'', ''y'' '
        df.index.name = 'cell_label'
        self._centroid = df.copy()

    # -------- METHODS -------- #
    def ini_centroids(self):
        d = {
            'x': self.ini_cell_props['x0'],
            'y': self.ini_cell_props['y0'],
        }
        df = pd.DataFrame(d)
        return df.copy()

    def ini_cov(self):
        mcr = self.dapi_mean_cell_radius()
        cov = mcr * mcr * np.eye(2, 2)
        return np.tile(cov, (self.nC, 1, 1))

    def dapi_mean_cell_radius(self):
        return np.nanmean(np.sqrt(self.ini_cell_props['area'] / np.pi)) * 0.5

    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.yx_coords)
        return nbrs

    def scatter_matrix(self, spots):
        mu_bar = self.centroid.values
        prob = spots.parent_cell_prob[:, :-1]
        id = spots.parent_cell_id[:, :-1]
        xy_spots = spots.xy_coords
        out = self.ini_cov() * self.nu_0

        mu_x = mu_bar[id, 0]  # array of size [nS, N] with the x-coord of the centroid of the N closest cells
        mu_y = mu_bar[id, 1]  # array of size [nS, N] with the y-coord of the centroid of the N closest cells

        N = mu_x.shape[1]
        _x = np.tile(xy_spots[:, 0], (N, 1)).T  # array of size [nS, N] populated with the x-coord of the spot
        _y = np.tile(xy_spots[:, 1], (N, 1)).T  # array of size [nS, N] populated with the y-coord of the spot
        x_centered = _x - mu_x  # subtract the cell centroid x-coord from the spot x-coord
        y_centered = _y - mu_y  # subtract the cell centroid y-coord from the spot x-coord

        el_00 = prob * x_centered * x_centered  # contribution to the scatter matrix's [0, 0] element
        el_11 = prob * y_centered * y_centered  # contribution to the scatter matrix's off-diagonal element
        el_01 = prob * x_centered * y_centered  # contribution to the scatter matrix's [1, 1] element

        # Aggregate all contributions to get the scatter matrix
        agg_00 = npg.aggregate(id.ravel(), el_00.ravel(), size=self.nC)
        agg_11 = npg.aggregate(id.ravel(), el_11.ravel(), size=self.nC)
        agg_01 = npg.aggregate(id.ravel(), el_01.ravel(), size=self.nC)

        # Return now the scatter matrix. Some cell might not have any spots nearby. For those empty cells,
        # the scatter matrix will be a squared zero array. That is fine.
        out[:, 0, 0] = agg_00
        out[:, 1, 1] = agg_11
        out[:, 0, 1] = agg_01
        out[:, 1, 0] = agg_01

        return out

    def read_image_objects(self, img_obj, cfg):
        meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
        relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

        # append 1 for the misreads
        relCellRadius = np.append(1, relCellRadius)

        nom = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(cfg['InsideCellBonus'])) + np.exp(cfg['InsideCellBonus'])
        denom = np.exp(-0.5) * (1 - np.exp(cfg['InsideCellBonus'])) + np.exp(cfg['InsideCellBonus'])
        CellAreaFactor = nom / denom

        out = {}
        out['area_factor'] = CellAreaFactor
        # out['area_factor'] = np.ones(CellAreaFactor.shape)
        # logger.info('Overriden CellAreaFactor = 1')
        out['rel_radius'] = relCellRadius
        out['area'] = np.append(np.nan, img_obj.area)
        out['x0'] = np.append(-sys.maxsize, img_obj.x0.values)
        out['y0'] = np.append(-sys.maxsize, img_obj.y0.values)
        out['cell_label'] = np.append(0, img_obj.label.values)
        if 'old_label' in img_obj.columns:
            out['cell_label_old'] = np.append(0, img_obj.old_label.values)
        # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
        # and will be used to get all the misreads. It was given the label=0 and some very small
        # negative coords

        return out, meanCellRadius

# ----------------------------------------Class: Genes--------------------------------------------------- #
class Genes(object):
    def __init__(self, spots):
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta_bar = None
        self._logeta_bar = None
        self.nG = len(self.gene_panel)

    @property
    def eta_bar(self):
        return self._eta_bar

    @property
    def logeta_bar(self):
        return self._logeta_bar

    def init_eta(self, a, b):
        self._eta_bar = np.ones(self.nG) * (a / b)
        self._logeta_bar = np.ones(self.nG) * self._digamma(a, b)

    def calc_eta(self, a, b):
        self._eta_bar = a/b
        self._logeta_bar = self._digamma(a, b)

    def _digamma(self, a, b):
        return scipy.special.psi(a) - np.log(b)


# ----------------------------------------Class: Spots--------------------------------------------------- #
class Spots(object):
    def __init__(self, spots_df, config):
        self._parent_cell_prob = None
        self._parent_cell_id = None
        self.config = config
        self.data = self.read(spots_df)
        self.nS = self.data.shape[0]
        self.unique_gene_names = None
        self._gamma_bar = None
        self._log_gamma_bar = None
        [_, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True, return_counts=True)

    def __getstate__(self):
        # set here attributes to be excluded from serialisation (pickling)
        # It makes the pickle filesize smaller but maybe this will have to
        # change in the future.
        # These two attributes take up a lot of space on the disk:
        # _gamma_bar and _log_gamma_bar
        # FYI: https://realpython.com/python-pickle-module/
        attributes = self.__dict__.copy()
        del attributes['_gamma_bar']
        del attributes['_log_gamma_bar']
        return attributes

    # -------- PROPERTIES -------- #
    @property
    def gamma_bar(self):
        return self._gamma_bar.astype(self.config['dtype'])

    @gamma_bar.setter
    def gamma_bar(self, val):
        self._gamma_bar = val.astype(self.config['dtype'])

    @property
    def log_gamma_bar(self):
        return self._log_gamma_bar

    @log_gamma_bar.setter
    def log_gamma_bar(self, val):
        self._log_gamma_bar = val

    @property
    def xy_coords(self):
        lst = list(zip(*[self.data.x, self.data.y]))
        return np.array(lst)

    @property
    def parent_cell_prob(self):
        return self._parent_cell_prob

    @parent_cell_prob.setter
    def parent_cell_prob(self, val):
        self._parent_cell_prob = val

    @property
    def parent_cell_id(self):
        return self._parent_cell_id

    @parent_cell_id.setter
    def parent_cell_id(self, val):
        self._parent_cell_id = val


    # -------- METHODS -------- #
    def read(self, spots_df):
        # No need for x_global, y_global to be in the spots_df at first place.
        # Instead of renaming here, you could just use 'x' and 'y' when you
        # created the spots_df
        spots_df = spots_df.rename(columns={'x_global': 'x', 'y_global': 'y'})

        # remove a gene if it is on the exclude list
        exclude_genes = self.config['exclude_genes']
        gene_mask = [True if d not in exclude_genes else False for d in spots_df.target]
        spots_df = spots_df.loc[gene_mask]
        return spots_df.rename_axis('spot_id').rename(columns={'target': 'gene_name'})

    def cells_nearby(self, cells: Cells) -> np.array:
        spotYX = self.data[['y', 'x']].values

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotYX)

        # last column is for misreads.
        neighbors[:, -1] = 0
        return neighbors

    def ini_cellProb(self, neighbors, cfg):
        nS = self.data.shape[0]
        nN = cfg['nNeighbors'] + 1
        SpotInCell = self.data.label
        # assert (np.all(SpotInCell.index == neighbors.index))

        ## sanity check (this actually needs to be rewritten)
        # mask = np.greater(SpotInCell, 0, where=~np.isnan(SpotInCell))
        # sanity_check = neighbors[mask, 0] + 1 == SpotInCell[mask]
        # assert ~any(sanity_check), "a spot is in a cell not closest neighbor!"

        pSpotNeighb = np.zeros([nS, nN])
        pSpotNeighb[neighbors == SpotInCell.values[:, None]] = 1
        pSpotNeighb[SpotInCell == 0, -1] = 1

        ## Add a couple of checks here
        return pSpotNeighb

    def loglik(self, cells, cfg):
        # area = cells.ini_cell_props['area'][1:]
        # mcr = np.mean(np.sqrt(area / np.pi)) * 0.5  # This is the meanCellRadius
        mcr = cells.mcr
        dim = 2  # dimensions of the normal distribution: Bivariate
        # Assume a bivariate normal and calc the likelihood
        D = -self.Dist ** 2 / (2 * mcr ** 2) - dim/2 * np.log(2 * np.pi * mcr ** 2)

        # last column (nN-closest) keeps the misreads,
        D[:, -1] = np.log(cfg['MisreadDensity'])

        mask = np.greater(self.data.label.values, 0, where=~np.isnan(self.data.label.values))
        D[mask, 0] = D[mask, 0] + cfg['InsideCellBonus']
        return D

    def zero_class_counts(self, geneNo, pCellZero):
        """
        Gene counts for the zero expressing class
        """
        # for each spot get the ids of the 3 nearest cells
        spotNeighbours = self.parent_cell_id[:, :-1]

        # get the corresponding probabilities
        neighbourProb = self.parent_cell_prob[:, :-1]

        # prob that a spot belongs to a zero expressing cell
        pSpotZero = np.sum(neighbourProb * pCellZero[spotNeighbours], axis=1)

        # aggregate per gene id
        TotPredictedZ = np.bincount(geneNo, pSpotZero)
        return TotPredictedZ

    def gammaExpectation(self, rho, beta):
        '''
        :param r:
        :param b:
        :return: Expectetation of a rv X following a Gamma(r,b) distribution with pdf
        f(x;\alpha ,\beta )= \frac{\beta^r}{\Gamma(r)} x^{r-1}e^{-\beta x}
        '''

        # sanity check
        # assert (np.all(rho.coords['cell_id'].data == beta.coords['cell_id'])), 'rho and beta are not aligned'
        # assert (np.all(rho.coords['gene_name'].data == beta.coords['gene_name'])), 'rho and beta are not aligned'

        dtype = self.config['dtype']
        r = rho[:, :, None]
        if dtype == np.float64:
            gamma = np.empty(beta.shape)
            ne.evaluate('r/beta', out=gamma)
            return gamma
        else:
            return (r/beta).astype(dtype)

    def logGammaExpectation(self, rho, beta):
        dtype = self.config['dtype']
        r = rho[:, :, None].astype(dtype)
        if dtype == np.float64:
            logb = np.empty(beta.shape)
            ne.evaluate("log(beta)", out=logb)
            return scipy.special.psi(r) - logb
        else:
            return scipy.special.psi(r) - np.log(beta).astype(dtype)


# ----------------------------------------Class: SingleCell--------------------------------------------------- #
class SingleCell(object):
    def __init__(self, scdata: pd.DataFrame, genes: np.array, config):
        self.isMissing = None  # Will be set to False is single cell data are assumed known and given as an input
                               # otherwise, if they are unknown, this will be set to True and the algorithm will
                               # try to estimate them
        # self.raw_data = self._raw_data(scdata, genes)
        self.config = config
        self._mean_expression, self._log_mean_expression = self._setup(scdata, genes, self.config)

    def _setup(self, scdata, genes, config):
        """
        calcs the mean (and the log-mean) gene counts per cell type. Note that
        some hyperparameter values have been applied before those means are derived.
        These hyperparameters and some bacic cleaning takes part in the functions
        called herein
        """
        if scdata is None:
            logger.info('Single Cell data are missing. Cannot determine meam expression per cell class.')
            logger.info('We will try to estimate the array instead')
            logger.info('Starting point is a diagonal array of size numGenes-by-numGenes')
            # expr = self._naive(scdata, genes)
            expr = self._diag(genes)
            self.isMissing = True
        else:
            expr = self._raw_data(scdata, genes)
            self.isMissing = False

        self.raw_data = expr
        me, lme = self._helper(expr.copy())
        dtype = self.config['dtype']

        assert me.columns[-1] == 'Zero', "Last column should be the Zero class"
        assert lme.columns[-1] == 'Zero', "Last column should be the Zero class"
        return me.astype(dtype), lme.astype(dtype)

    # -------- PROPERTIES -------- #
    @property
    def mean_expression(self):
        assert self._mean_expression.columns[-1] == 'Zero', "Last column should be the Zero class"
        return self._mean_expression

    @property
    def log_mean_expression(self):
        assert self._log_mean_expression.columns[-1] == 'Zero', "Last column should be the Zero class"
        return self._log_mean_expression

    @property
    def genes(self):
        return self.mean_expression.index.values

    @property
    def classes(self):
        return self.mean_expression.columns.values

    ## Helper functions ##
    def _set_axes(self, df):
        df = df.rename_axis("class_name", axis="columns").rename_axis('gene_name')
        return df

    def _remove_zero_cols(self, df):
        """
        Removes zero columns (ie if a column is populated by zeros only, then it is removed)
        :param da:
        :return:
        """
        out = df.loc[:, (df != 0).any(axis=0)]
        return out

    def _helper(self, expr):
        # order by column name
        # expr = expr.copy().sort_index(axis=0).sort_index(axis=1, key=lambda x: x.str.lower())
        # expr = expr.copy().sort_index(axis=0).sort_index(axis=1, key=natsort_keygen())
        expr = expr.copy().sort_index(axis=0).sort_index(axis=1, key=natsort_keygen(key=lambda y: y.str.lower()))

        # append at the end the Zero class
        expr['Zero'] = np.zeros([expr.shape[0], 1])

        # expr = self.config['Inefficiency'] * arr
        me = expr.rename_axis('gene_name').rename_axis("class_name", axis="columns")  # mean expression
        me = me + self.config['SpotReg']  # add the regularization parameter
        lme = np.log(me)  # log mean expression
        return me, lme

    def _gene_expressions(self, fitted, scale):
        """
        Finds the expected mean gene counts. The prior *IS NOT* taken
        into account. We use data evidence only
        For the zero class only the prior is used *AND NOT* data
        evidence.
        """

        # the prior on mean expression follows a Gamma(m * M , m), where M is the starting point (the initial
        # array) of single cell data
        # 07-May-2023. Hiding m from the config.py. Should bring it back at a later version
        # m = self.config['m']
        m = 1

        a = fitted + m * (self.raw_data + self.config['SpotReg'])
        b = scale + m
        me = a / b
        lme = scipy.special.psi(a) - np.log(b)

        # the expressions for the zero class are a 0.0 plus the regularition param
        zero_col = np.zeros(me.shape[0]) + self.config['SpotReg']
        me = me.assign(Zero=zero_col)
        # For the mean of the log-expressions, again only the prior is used for the Zero class
        zero_col_2 = scipy.special.psi(m * zero_col) - np.log(m)
        lme = lme.assign(Zero=zero_col_2)
        return me, lme

    def _keep_labels_unique(self, scdata):
        """
        In the single cell data you might find cases where two or more rows have the same gene label
        In these cases keep the row with the highest total gene count
        """

        # 1. get the row total and assign it to a new column
        scdata = scdata.assign(total=scdata.sum(axis=1))

        # 2. rank by gene label and total gene count and keep the one with the highest total
        scdata = scdata.sort_values(['gene_name', 'total'], ascending=[True, False]).groupby('gene_name').head(1)

        # 3. Drop the total column and return
        return scdata.drop(['total'], axis=1)

    def _raw_data(self, scdata, genes):
        """
        Reads the raw single data, filters out any genes outside the gene panel and then it
        groups by the cell type
        """
        assert np.all(scdata >= 0), "Single cell dataframe has negative values"
        logger.info(' Single cell data passed-in have %d genes and %d cells' % (scdata.shape[0], scdata.shape[1]))

        logger.info(' Single cell data: Keeping counts for the gene panel of %d only' % len(genes))
        df = scdata.loc[genes]

        # set the axes labels
        df = self._set_axes(df)

        # remove any rows with the same gene label
        df = self._keep_labels_unique(df)

        df = self._remove_zero_cols(df.copy())
        dfT = df.T

        logger.info(' Single cell data: Grouping gene counts by cell type. Aggregating function is the mean.')
        out = dfT.groupby(dfT.index.values).agg('mean').T
        logger.info(' Grouped single cell data have %d genes and %d cell types' % (out.shape[0], out.shape[1]))
        return out

    def _diag(self, genes):
        # logger.info('******************************************************')
        # logger.info('*************** DIAGONAL SINGLE CELL DATA ***************')
        # logger.info('******************************************************')
        nG = len(genes)
        mgc = self.config['mean_gene_counts_per_class']  # the avg gene count per cell. Better expose that so it can be set by the user.
        arr = mgc * np.eye(nG)
        labels = ['class_%d' % (i+1) for i, _ in enumerate(genes)]
        df = pd.DataFrame(arr).set_index(genes)
        df.columns = labels
        return df

# ---------------------------------------- Class: CellType --------------------------------------------------- #
class CellType(object):
    def __init__(self, single_cell):
        assert single_cell.classes[-1] == 'Zero', "Last label should be the Zero class"
        self._names = single_cell.classes
        self._alpha = None
        self.single_cell_data_missing = single_cell.isMissing

    @property
    def names(self):
        assert self._names[-1] == 'Zero', "Last label should be the Zero class"
        return self._names

    @property
    def nK(self):
        return len(self.names)

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        self._alpha = val

    @property
    def pi_bar(self):
        return self.alpha / self.alpha.sum()

    @property
    def logpi_bar(self):
        return scipy.special.psi(self.alpha) - scipy.special.psi(self.alpha.sum())

    @property
    def prior(self):
        return self.pi_bar

    @property
    def log_prior(self):
        if self.single_cell_data_missing:
            return self.logpi_bar
        else:
            return np.log(self.prior)

    def size(self, cells):
        """
        calcs the size of a cell class, ie how many members (ie cells) each cell type has
        """
        return cells.classProb.sum(axis=0)

    def ini_prior(self):
        self.alpha = self.ini_alpha()

    def ini_alpha(self):
        return np.append(np.ones(self.nK - 1), sum(np.ones(self.nK - 1)))





