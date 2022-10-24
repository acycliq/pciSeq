import sys
import numpy as np
import pandas as pd
import datetime
import numpy_groupies as npg
import scipy
import numexpr as ne
from sklearn.neighbors import NearestNeighbors
from pciSeq.src.cell_call.log_config import logger


class Cells(object):
    # Get rid of the properties where not necessary!!
    def __init__(self, _cells_df, config):
        self.config = config
        self.cell_props, self.mcr = self.read_image_objects(_cells_df, config)
        self.nC = len(self.cell_props['cell_label'])
        self.classProb = None
        self.class_names = None
        self._prior = None
        self._cov = self.ini_cov()
        self.nu_0 = 30  # need to move that into config.py. Degrees of freedom of the Wishart prior. Use the mean gene counts per cell to set nu_0
        # self.rho_1 = self.config['relax_segmentation: rho_1']
        # self.rho_2 = self.config['relax_segmentation: rho_2']
        self._centroid = self.ini_centroids()
        self._gene_counts = None
        self._background_counts = None
        self._alpha = None

    # -------- PROPERTIES -------- #
    @property
    def zyx_coords(self):
        coords = [d for d in zip(self.cell_props['z'], self.cell_props['y'], self.cell_props['x']) if not np.isnan(d).any()]
        return np.array(coords)

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
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        self._alpha = val

    @property
    def centroid(self):
        # lst = list(zip(*[self._centroid['x'], self._centroid['y']]))
        return self._centroid.copy()

    @centroid.setter
    def centroid(self, df):
        assert isinstance(df, pd.DataFrame), 'Input should be a dataframe'
        assert set(df.columns.values) == {'x', 'y', 'z'}, 'Dataframe columns should be ''x'', ''y'' and ''z'' '
        self._centroid = df.copy()

    @property
    def cov(self):
        return self._cov

    @cov.setter
    def cov(self, val):
        self._cov = val

    @property
    def sigma_x(self):
        return np.sqrt(self.cov[:, 0, 0])

    @property
    def sigma_y(self):
        return np.sqrt(self.cov[:, 1, 1])

    @property
    def sigma_z(self):
        return np.sqrt(self.cov[:, 2, 2])

    @property
    def corr(self):
        sigma_x = np.sqrt(self.cov[:, 0, 0])
        sigma_y = np.sqrt(self.cov[:, 1, 1])
        sigma_z = np.sqrt(self.cov[:, 2, 2])
        cov_xy = self.cov[:, 0, 1]
        cov_xz = self.cov[:, 0, 2]
        cov_yz = self.cov[:, 1, 2]
        rho_xy = cov_xy / (sigma_x * sigma_y)
        rho_xz = cov_xz / (sigma_x * sigma_z)
        rho_yz = cov_yz / (sigma_y * sigma_z)
        return [rho_xy, rho_xz, rho_yz]

    @property
    def ellipsoid_attributes(self):
        """
        convenience property that returns a list with
        the ellipoid's centroid, correlation and standard
        deviation across the x-axis and y-axis
        """
        mu = self.centroid.values.tolist()
        sigma_x = self.sigma_x.tolist()
        sigma_y = self.sigma_y.tolist()
        rho = self.corr
        return list(zip(mu, rho, sigma_x, sigma_y))

    # -------- METHODS -------- #
    def ini_centroids(self):
        d = {
            'x': self.cell_props['x'],
            'y': self.cell_props['y'],
            'z': self.cell_props['z'],
        }
        df = pd.DataFrame(d)
        return df.copy()

    def ini_cov(self):
        dim = 3
        cov = self.mcr * self.mcr * np.eye(dim, dim)
        return np.tile(cov, (self.nC, 1, 1))

    def dapi_mean_cell_radius(self):
        return np.nanmean(np.sqrt(self.cell_props['area'] / np.pi)) * 0.5

    def nn(self):
        n = self.config['nNeighbors'] + 1
        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = NearestNeighbors(n_neighbors=n, algorithm='ball_tree').fit(self.zyx_coords)
        return nbrs

    # def geneCountsPerKlass(self, single_cell_data, egamma, ini):
    #     # ********************************************
    #     # DEPRECATED. Replaced by a simple einsum call
    #     # ********************************************
    #     temp = np.einsum('ck, c, cgk -> gk', self.classProb, self.cell_props['area_factor'], egamma)
    #
    #     # total counts predicted by all cells of each class (nG, nK)
    #     ClassTotPredicted = temp * (single_cell_data.mean_expression + ini['SpotReg'])
    #
    #     # total of each gene
    #     isZero = ClassTotPredicted.columns == 'Zero'
    #     labels = ClassTotPredicted.columns.values[~isZero]
    #     TotPredicted = ClassTotPredicted[labels].sum(axis=1)
    #     return TotPredicted

    def scatter_matrix(self, spots):
        mu_bar = self.centroid.values
        prob = spots.parent_cell_prob[:, :-1]
        id = spots.parent_cell_id[:, :-1]
        xyz_spots = spots.xyz_coords
        out = self.ini_cov() * self.nu_0

        mu_x = mu_bar[id, 0]  # array of size [nS, N] with the x-coord of the centroid of the N closest cells
        mu_y = mu_bar[id, 1]  # array of size [nS, N] with the y-coord of the centroid of the N closest cells
        mu_z = mu_bar[id, 2]  # array of size [nS, N] with the z-coord of the centroid of the N closest cells

        N = mu_x.shape[1]
        _x = np.tile(xyz_spots[:, 0], (N, 1)).T  # array of size [nS, N] populated with the x-coord of the spot
        _y = np.tile(xyz_spots[:, 1], (N, 1)).T  # array of size [nS, N] populated with the y-coord of the spot
        _z = np.tile(xyz_spots[:, 2], (N, 1)).T  # array of size [nS, N] populated with the z-coord of the spot

        x_centered = _x - mu_x  # subtract the cell centroid x-coord from the spot x-coord
        y_centered = _y - mu_y  # subtract the cell centroid y-coord from the spot y-coord
        z_centered = _z - mu_z  # subtract the cell centroid z-coord from the spot z-coord

        el_00 = prob * x_centered * x_centered  # contribution to the scatter matrix's [0, 0] element
        el_11 = prob * y_centered * y_centered  # contribution to the scatter matrix's [1, 1] element
        el_22 = prob * z_centered * z_centered  # contribution to the scatter matrix's [2, 2] element

        el_01 = prob * x_centered * y_centered  # contribution to the scatter matrix's [0, 1] element
        el_02 = prob * x_centered * z_centered  # contribution to the scatter matrix's [0, 2] element
        el_12 = prob * y_centered * z_centered  # contribution to the scatter matrix's [1, 2] element

        # Aggregate all contributions to get the scatter matrix
        agg_00 = npg.aggregate(id.ravel(), el_00.ravel(), size=self.nC)
        agg_11 = npg.aggregate(id.ravel(), el_11.ravel(), size=self.nC)
        agg_22 = npg.aggregate(id.ravel(), el_22.ravel(), size=self.nC)

        agg_01 = npg.aggregate(id.ravel(), el_01.ravel(), size=self.nC)
        agg_02 = npg.aggregate(id.ravel(), el_02.ravel(), size=self.nC)
        agg_12 = npg.aggregate(id.ravel(), el_12.ravel(), size=self.nC)

        # Return now the scatter matrix. Some cell might not have any spots nearby. For those empty cells,
        # the scatter matrix will be a squared zero array. That is fine.
        out[:, 0, 0] = agg_00
        out[:, 1, 1] = agg_11
        out[:, 2, 2] = agg_22

        out[:, 0, 1] = agg_01
        out[:, 0, 2] = agg_02
        out[:, 1, 2] = agg_12

        out[:, 1, 0] = agg_01
        out[:, 2, 0] = agg_02
        out[:, 2, 1] = agg_12

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
        if cfg['relax_segmentation'] or cfg['is_3D']:
            out['area_factor'] = np.ones(CellAreaFactor.shape[0])
        else:
            out['area_factor'] = CellAreaFactor
        # out['area_factor'] = np.ones(CellAreaFactor.shape)
        # logger.info('Overriden CellAreaFactor = 1')
        out['rel_radius'] = relCellRadius
        out['area'] = np.append(np.nan, img_obj.area)
        out['x'] = np.append(-sys.maxsize, img_obj.x.values)
        out['y'] = np.append(-sys.maxsize, img_obj.y.values)
        out['z'] = np.append(-sys.maxsize, img_obj.z.values)
        out['cell_label'] = np.append(0, img_obj.label.values)
        if 'old_label' in img_obj.columns:
            out['cell_label_old'] = np.append(0, img_obj.old_label.values)
        # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
        # and will be used to get all the misreads. It was given the label=0 and some very small
        # negative coords

        return out, meanCellRadius


    def export_class_prob(self, conn, iter_num, has_converged):
        pass



# ----------------------------------------Class: Genes--------------------------------------------------- #
class Genes(object):
    def __init__(self, spots):
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta_bar = None
        self._logeta_bar = None
        self.nG = len(self.gene_panel)
        # self.con = con
        # self.con.executescript("CREATE TABLE gene_efficiency (gene varchar,"
        #                        "gene_efficiency,"
        #                        "iteration,"
        #                        "has_converged,utc"
        #                        ")"
        #                        )
        # self.con.executescript('''
        #     CREATE TRIGGER gene_efficiency_row
        #     AFTER INSERT ON gene_efficiency
        #     WHEN (SELECT count(*) FROM gene_efficiency) > 0
        #     BEGIN
        #         SELECT RAISE (FAIL, 'full');
        #     END;
        # ''')
        # self.con.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_gene_iteration ON gene_efficiency("gene", "iteration");')

    @property
    def eta(self):
        raise Exception

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

    def db_save(self, con, run, converged, db_opts):
        self.db_gene_efficiency(con, run, converged, db_opts)

    def db_gene_efficiency(self, con, run, converged, db_opts):
        df = pd.DataFrame(data=self.eta_bar, index=self.gene_panel, columns=['gene_efficiency'])
        df.index.name = 'gene'
        df['iteration'] = run
        df['has_converged'] = converged
        df['utc'] = datetime.datetime.utcnow()
        df = df.reset_index()
        df.to_sql(name='gene_efficiency', con=con, if_exists=db_opts['if_table_exists'], index=False)
        con.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_gene_iteration ON gene_efficiency("gene", "iteration");')


# ----------------------------------------Class: Spots--------------------------------------------------- #
class Spots(object):
    def __init__(self, spots_df, config):
        self._parent_cell_prob = None
        self._parent_cell_id = None
        self.config = config
        self.data, self.data_excluded = self.read(spots_df)
        self.nS = self.data.shape[0]
        self.call = None
        self.unique_gene_names = None
        self._gamma_bar = None
        self._log_gamma_bar = None
        [_, self.gene_id, self.counts_per_gene] = np.unique(self.data.gene_name.values, return_inverse=True, return_counts=True)

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
    def xyz_coords(self):
        lst = list(zip(*[self.data.x, self.data.y, self.data.z]))
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
        # spots_df = spots_df.rename(columns={'x_global': 'x', 'y_global': 'y', 'z_global': 'z'})

        # remove a gene if it is on the exclude list
        exclude_genes = self.config['exclude_genes']
        gene_mask = [True if d not in exclude_genes else False for d in spots_df.gene_name]
        neg_gene_mask = [True if d in exclude_genes else False for d in spots_df.gene_name]
        spots_copy = spots_df.copy()
        spots_df = spots_copy.loc[gene_mask]
        spots_excluded_df = spots_copy.loc[neg_gene_mask]
        return spots_df, spots_excluded_df

    def cells_nearby(self, cells: Cells) -> np.array:
        spotZYX = self.data[['z', 'y', 'x']]

        # for each spot find the closest cell (in fact the top nN-closest cells...)
        nbrs = cells.nn()
        self.Dist, neighbors = nbrs.kneighbors(spotZYX.values)

        # last column is for misreads.
        neighbors[:, -1] = 0

        # make an array assigning 100% prob of any given cell belonging to its closest neighbour
        cellProb = np.zeros(neighbors.shape, dtype=np.uint32)
        cellProb[:, 0] = np.ones(neighbors.shape[0])
        return neighbors, cellProb

    def ini_cellProb(self, neighbors, cfg):
        nS = self.data.shape[0]
        nN = cfg['nNeighbors'] + 1
        SpotInCell = self.data.label
        # assert (np.all(SpotInCell.index == neighbors.index))

        # sanity check (this actually needs to be rewritten)
        mask = np.greater(SpotInCell, 0, where=~np.isnan(SpotInCell))
        sanity_check = neighbors[mask, 0] + 1 == SpotInCell[mask]
        assert ~any(sanity_check), "a spot is in a cell not closest neighbor!"

        pSpotNeighb = np.zeros([nS, nN])
        pSpotNeighb[neighbors == SpotInCell.values[:, None]] = 1
        pSpotNeighb[SpotInCell == 0, -1] = 1

        ## Add a couple of checks here
        return pSpotNeighb

    def loglik(self, cells, cfg):
        # area = cells.cell_props['area'][1:]
        # mcr = np.mean(np.sqrt(area / np.pi)) * 0.5  # This is the meanCellRadius
        mcr = cells.mcr
        dim = 2  # dimensions of the normal distribution: Bivariate
        # Assume a bivariate normal and calc the likelihood
        D = -self.Dist ** 2 / (2 * mcr ** 2) - dim/2 * np.log(2 * np.pi * mcr ** 2)

        # last column (nN-closest) keeps the misreads,
        D[:, -1] = np.log(cfg['MisreadDensity'])

        mask = np.greater(self.data.label, 0, where=~np.isnan(self.data.label))
        D[mask, 0] = D[mask, 0] + cfg['InsideCellBonus']
        return D

    def mvn_loglik(self, data, cell_label, cells):
        centroids = cells.centroid.values[cell_label]
        covs = cells.cov[cell_label]
        param = list(zip(*[data, centroids, covs]))
        # out = [self.loglik_contr(p) for i, p in enumerate(param)]
        # out_2 = [multivariate_normal.logpdf(p[0], p[1], p[2]) for i, p in enumerate(param)]
        out = self.multiple_logpdfs(data, centroids, covs)
        return out

    def multiple_logpdfs(self, x, means, covs):
        """
        vectorised mvn log likelihood evaluated at multiple pairs of (centroid_1, cov_1), ..., (centroid_N, cov_N)
        Taken from http://gregorygundersen.com/blog/2020/12/12/group-multivariate-normal-pdf/
        """
        # Thankfully, NumPy broadcasts `eigh`.
        vals, vecs = np.linalg.eigh(covs)

        # Compute the log determinants across the second axis.
        logdets = np.sum(np.log(vals), axis=1)

        # Invert the eigenvalues.
        valsinvs = 1. / vals

        # Add a dimension to `valsinvs` so that NumPy broadcasts appropriately.
        Us = vecs * np.sqrt(valsinvs)[:, None]
        devs = x - means

        # Use `einsum` for matrix-vector multiplications across the first dimension.
        devUs = np.einsum('ni,nij->nj', devs, Us)

        # Compute the Mahalanobis distance by squaring each term and summing.
        mahas = np.sum(np.square(devUs), axis=1)

        # Compute and broadcast scalar normalizers.
        dim = len(vals[0])
        log2pi = np.log(2 * np.pi)

        return -0.5 * (dim * log2pi + mahas + logdets)

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

    def db_save(self, conn):
        self.db_save_spots(conn)

    def db_save_spots(self, con):
        # spots are persistent, they do not change. Just save it once
        try:
            df = pd.DataFrame(data=self.data[['x', 'y', 'z', 'label', 'gene_name']], columns=['x', 'y', 'z', 'label', 'gene_name'])
            df = df.set_index('gene_name')
            df['utc'] = datetime.datetime.utcnow()
            df = df.reset_index()
            df.to_sql(name='spots', con=con, if_exists='fail', index=False)
            con.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_class_iteration ON spots("gene_name");')
        except ValueError:
            pass


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
            logger.info('Single Cell data are missing. We will try to estimate it')
            logger.info('Starting point is a diagonal array')
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
        expr = expr.copy().sort_index(axis=0).sort_index(axis=1, key=lambda x: x.str.lower())

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
        Ffor the zero class only the prior is used *AND NOT* data
        evidence.
        """
        m = self.config['m']
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
        logger.info('******************************************************')
        logger.info('*************** DIAGONAL SINGLE CELL DATA ***************')
        logger.info('******************************************************')
        nG = len(genes)
        # nK = 71
        mgc = 15
        arr = mgc * np.eye(nG)
        # labels = ["label_" + str(d) for d in range(arr.shape[1])]
        labels = [s + '_class' for s in genes]
        df = pd.DataFrame(arr).set_index(genes)
        df.columns = labels

        return df

    def db_save(self, conn, run, converged, db_opts):
        self.db_save_mean_expressions(conn, run, converged, db_opts)

    def db_save_mean_expressions(self, con, run, converged, db_opts):
        """
        pushing the average expressions to the db. Most of the times, these are persistent data, therefore
        appending the columns, 'iteration' and 'has_converged' will not make a lot of sense unless single
        cell data are missing, hence they will be estimated on the fly
        """
        df = self.mean_expression.copy()
        df['utc'] = datetime.datetime.utcnow()
        df['iteration'] = run
        df['has_converged'] = converged
        df = df.reset_index()
        df.to_sql(name='mean_expression', con=con, if_exists=db_opts['if_table_exists'], index=False)
        con.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_gene_iteration ON mean_expression("gene_name", "iteration");')



# ---------------------------------------- Class: CellType --------------------------------------------------- #
class CellType(object):
    def __init__(self, single_cell, config):
        assert single_cell.classes[-1] == 'Zero', "Last label should be the Zero class"
        self._names = single_cell.classes
        self.config = config
        # self._prior = None
        self._alpha = None

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
        if not self.config['is_3D'] or not self.config['relax_segmentation']:
            assert set(self.pi_bar) == {0.5, 0.5 / (self.nK - 1)}
        return self.pi_bar

    @property
    def log_prior(self):
        if self.config['is_3D'] or self.config['relax_segmentation']:
            return self.logpi_bar
        else:
            assert set(self.pi_bar) == {0.5, 0.5 / (self.nK - 1)}
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

    def db_save(self, con, run, converged, db_opts):
        self.db_cell_type_prior(con, run, converged, db_opts)

    def db_cell_type_prior(self, con, run, converged, db_opts):
        df = pd.DataFrame(data=self.pi_bar, index=self.names, columns=['weight'])
        df.index.name = 'class'
        df['iteration'] = run
        df['has_converged'] = converged
        df['utc'] = datetime.datetime.utcnow()
        df = df.reset_index()
        df.to_sql(name='cell_type_prior', con=con, if_exists=db_opts['if_table_exists'], index=False)
        con.execute('CREATE UNIQUE INDEX IF NOT EXISTS ix_class_iteration ON cell_type_prior("class", "iteration");')




