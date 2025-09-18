# Standard library imports
import logging
from typing import Tuple, Dict

# Third party imports
import numpy as np
import pandas as pd
import xarray as xr
import scipy
from natsort import natsort_keygen

# Local imports
from ..utils.cell_utils import read_image_objects, keep_labels_unique
from ..utils.ops_utils import gene_density

singleCell_logger = logging.getLogger(__name__)


class SingleCell(object):
    """
    Handles single-cell RNA sequencing reference data, including mean
    expression levels per cell type. Supports integration of single-cell
    data into broader analyses.

    Attributes:
        isMissing (bool): Indicates if single-cell data is missing.
        config (dict): Configuration parameters for single-cell data.
        _mean_expression (pd.DataFrame): Mean expression levels.
        _log_mean_expression (pd.DataFrame): Log mean expression levels.
    """

    def __init__(self, scdata: pd.DataFrame, genes: np.ndarray, spots, config: Dict):
        """
        Initializes the SingleCell object with single-cell data and configuration.

        Parameters:
            scdata (pd.DataFrame): Single-cell data.
            genes (np.array): Array of gene names.
            config (dict): Configuration parameters for single-cell data.
        """
        self.isMissing = None  # Will be set to False if single cell data are assumed known and given as an input
        # otherwise, if they are unknown, this will be set to True and the algorithm will
        # try to estimate them
        # self.raw_data = self._raw_data(scdata, genes)
        self.config = config
        self.spots = spots
        self._mean_expression_adj, self._log_mean_expression_adj, self._mean_expression = self._setup(scdata, genes, self.config)


    def _setup(self, scdata: pd.DataFrame, genes: np.ndarray, config: Dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Sets up the mean and log mean expression levels.

        Parameters:
            scdata (pd.DataFrame): Single-cell data.
            genes (np.array): Array of gene names.
            config (dict): Configuration parameters.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Mean and log mean expression levels.
        """
        if scdata is None:
            singleCell_logger.info('Single Cell data are missing. Cannot determine mean expression per cell class.')
            singleCell_logger.info('We will try to estimate the array instead')
            singleCell_logger.info('Starting point is a diagonal array of size numGenes-by-numGenes')
            expr = self._diag(genes)
            self.isMissing = True
        else:
            expr = self._raw_data(scdata, genes)
            self.isMissing = False

        self.raw_data = expr

        # get the mean (and log-mean) expression data by cell type.
        # Figures have been scaled by the gene inefficiency
        me_adj, lme_adj, mean_expr = self._helper(expr.copy())

        assert me_adj.class_name.values[-1] == 'Zero', "Last column should be the Zero class"
        assert lme_adj.class_name.values[-1] == 'Zero', "Last column should be the Zero class"
        return me_adj.astype(np.float32), lme_adj.astype(np.float32), mean_expr.astype(np.float32)

    # -------- PROPERTIES -------- #
    @property
    def mean_expression_adj(self):
        """Returns the mean expression levels adjusted by the initial gene inefficiency."""
        assert self._mean_expression_adj.class_name.values[-1] == 'Zero', "Last column should be the Zero class"
        return self._mean_expression_adj

    @property
    def log_mean_expression(self):
        """Returns the log mean expression levels adjusted by the initial gene inefficiency."""
        assert self._log_mean_expression_adj.class_name.values[-1] == 'Zero', "Last column should be the Zero class"
        return self._log_mean_expression_adj

    @property
    def mean_expression(self):
        """Returns the mean gene counts per cell class"""
        return self._mean_expression

    @property
    def genes(self):
        """
        Returns the gene names.
        WARNING: you can get the gene names from the Genes object and its gene_panel property.
        They should be the same but is there any value having two places for the same thing?
        """
        return self.mean_expression_adj.index.values

    @property
    def classes(self):
        """Returns the class names."""
        return self.mean_expression_adj.class_name.values

    # Helper functions #
    def _set_axes(self, df):
        """
        Sets the axes labels for a DataFrame.
        """
        df = df.rename_axis("class_name", axis="columns").rename_axis('gene_name')
        return df

    def _remove_zero_cols(self, df):
        """
        Removes columns with all zero values from a DataFrame.
        """
        out = df.loc[:, (df != 0).any(axis=0)]
        return out

    def _helper(self, expr_in):
        """
        Helper function to process expression data.

        Parameters:
            expr (pd.DataFrame): Expression data.

        Returns:
            tuple: Processed mean and log-mean expression data.
        """

        # order by column name
        expr = expr_in.copy().sort_index(axis=0).sort_index(axis=1, key=natsort_keygen(key=lambda y: y.str.lower()))

        # append at the end the Zero class
        expr['Zero'] = np.zeros([expr.shape[0], 1])
        expr = expr.rename_axis('gene_name').rename_axis("class_name", axis="columns")

        # apply the inefficiency
        me = expr * self.config['Inefficiency']

        # plane adjustment. This is the gene density in each plane.
        # It will be used to scale the single cell data depending on the location of the cell centroid
        plane_adj = gene_density(self.spots, self.config)

        # print("***** REMOVE THIS REMOVE THIS REMOVE THIS REMOVE THIS *****")
        # plane_adj = np.ones(plane_adj.shape) # REMOVE THIS REMOVE THIS REMOVE THIS REMOVE THIS

        # if an element is zero, set it to 1. Effectively that means that if the gene is not present on that plane,
        # then dont make any adjustment to the single cell data
        plane_adj[plane_adj == 0] = 1

        me = np.einsum('gk,pg->pgk', me, plane_adj)

        me = xr.DataArray(
            me,
            dims=['plane', 'gene_name', 'class_name'],
            coords={
                'plane': np.arange(len(me)),
                'gene_name': expr.index.values,
                'class_name': expr.columns.values
            }
        )


        # log mean expression
        lme = np.log(me + self.config['SpotReg'])
        return me, lme, expr

    def _gene_expressions(self, fitted, scale):
        """
        Calculates expected mean gene counts. The prior *IS NOT* taken
        into account. We use data evidence only
        For the zero class only the prior is used *AND NOT* data
        evidence.

        Parameters:
            fitted (np.array): Fitted values.
            scale (np.array): Scale values.

        Returns:
            tuple: Mean and log-mean gene expressions.
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

    def _raw_data(self, scdata: pd.DataFrame, genes: np.ndarray) -> pd.DataFrame:
        """
        Processes raw single-cell data, filtering out any genes outside the gene panel and grouping by cell type.

        Parameters:
            scdata (pd.DataFrame): Single-cell data.
            genes (np.array): Array of gene names.

        Returns:
            pd.DataFrame: Processed single-cell data.
        """
        assert np.all(scdata >= 0), "Single cell dataframe has negative values"
        singleCell_logger.info('Single cell data passed-in have %d genes and %d cells' % (scdata.shape[0], scdata.shape[1]))
        singleCell_logger.info('Single cell data: Keeping counts for the gene panel of %d only' % len(genes))
        df = scdata.loc[genes]

        # set the axes labels
        df = self._set_axes(df)

        # remove any rows with the same gene label
        df = keep_labels_unique(df)

        df = self._remove_zero_cols(df.copy())
        dfT = df.T

        singleCell_logger.info('Single cell data: Grouping gene counts by cell type. Aggregating function is the mean.')
        out = dfT.groupby(dfT.index.values).agg('mean').T
        singleCell_logger.info('Grouped single cell data have %d genes and %d cell types' % (out.shape[0], out.shape[1]))
        return out

    def _diag(self, genes):
        """
        Creates a diagonal matrix for single-cell data initialization.

        Parameters:
            genes (np.array): Array of gene names.

        Returns:
            pd.DataFrame: Diagonal single-cell data.
        """
        nG = len(genes)
        mgc = self.config['mean_gene_counts_per_class']
        arr = mgc * np.eye(nG)
        labels = ['class_%d' % (i + 1) for i, _ in enumerate(genes)]
        df = pd.DataFrame(arr).set_index(genes)
        df.columns = labels
        return df

