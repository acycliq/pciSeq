# Standard library imports
import logging
from typing import Tuple, Dict, Any

# Third party imports
import numpy as np
import pandas as pd
import scipy

genes_logger = logging.getLogger(__name__)


class Genes(object):
    """
    Manages gene-specific data and calculations, including initialization and
    computation of gene expression parameters.

    Attributes:
        gene_panel (np.array): Array of unique gene names.
        _eta_bar (np.array): Eta bar values: This is basically the expected Gene inefficiency.
        _logeta_bar (np.array): Log eta bar values for genes.
        nG (int): Number of genes.
    """

    def __init__(self, spots, config: Dict):
        """
        Initializes the Genes object with spot data.

        Parameters:
            spots (Spots): Spots object containing spot data.
        """
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta_bar = None
        self._logeta_bar = None
        self.nG = len(self.gene_panel)
        self.config = config

    @property
    def eta_bar(self):
        """Returns the eta bar values for genes."""
        return self._eta_bar

    @property
    def logeta_bar(self):
        """Returns the log eta bar for genes (estimated mean of the posterior)."""
        return self._logeta_bar

    @property
    def inefficiency(self):
        """
        Returns the gene inefficiency
        The actual gene inefficiency is the estimated mean of the posterior (eta_bar)
        multiplied by the inefficiency (user-defined) value that was passed in the algo
        via the configuration file
        """
        return self.eta_bar * self.config['Inefficiency']

    def init_eta(self, a, b):
        """
        Initializes eta values for genes.

        Parameters:
            a (float): Parameter a for eta calculation.
            b (float): Parameter b for eta calculation.
        """
        self._eta_bar = np.ones(self.nG, dtype=np.float32) * (a / b)
        self._logeta_bar = np.ones(self.nG, dtype=np.float32) * self._digamma(a, b)

    def calc_eta(self, a, b):
        """
        Calculates eta values for genes.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.
        """
        a = a.astype(np.float32)
        b = b.astype(np.float32)
        self._eta_bar = a / b
        self._logeta_bar = self._digamma(a, b)

    def _digamma(self, a, b):
        """
        Calculates the digamma function for eta calculation.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.

        Returns:
            np.array: Digamma values.
        """
        return scipy.special.psi(a) - np.log(b)

    def get_inefficiency(self, gene=None):
        """
        Retrieve the inefficiency values for one or more genes.

        This is a convenience method that returns inefficiency values from the gene panel.
        it returns a  DataFrame with the genes and the corresponding inefficiencies.

        Parameters:
        ----------
        gene : str, list of str, or None (default: None)
            - If None, returns inefficiency values for all genes.
            - If a string, returns inefficiency for the specified gene as a DataFrame.
            - If a list of strings, returns inefficiency values for the specified genes.

        Returns:
        -------
        pandas.DataFrame
            A DataFrame with genes as the index and inefficiency values as the column.
            Missing genes will be included with NaN values.

        Raises:
        ------
        TypeError
            If `gene` is not a string, list of strings, or None.
        """
        df = pd.DataFrame(
            {'inefficiency': self.inefficiency},
            index=self.gene_panel
        )

        if gene is None:
            return df  # Return full DataFrame

        if isinstance(gene, str):
            return df.loc[[gene]] if gene in df.index else pd.DataFrame(columns=df.columns, index=[gene])

        if isinstance(gene, list):
            return df.reindex(gene)  # Handles missing genes gracefully (NaN for missing ones)

        raise TypeError("Expected gene to be a string, list, or None.")

