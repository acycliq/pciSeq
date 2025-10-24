# Standard library imports
import logging
from typing import Dict

# Third party imports
import numpy as np
import scipy

from .singleCell import SingleCell
from .cells import Cells

cellType_logger = logging.getLogger(__name__)


class CellClass(object):
    """
    Manages cell type classification, including prior probabilities and
    class assignments. Helps in understanding the distribution of different
    cell types and their characteristics within a dataset.

    Attributes:
        _names (np.array): Names of cell types.
        _alpha (np.array): Alpha values for cell types.
        config (dict): Configuration parameters for cell types.
        single_cell_data_missing (bool): Indicates if single-cell data is missing.
    """

    def __init__(self, single_cell: SingleCell, config: Dict):
        """
        Initializes the CellType object with single-cell data and configuration.

        Parameters:
            single_cell (SingleCell): SingleCell object containing single-cell data.
            config (dict): Configuration parameters for cell types.
        """
        assert single_cell.classes[-1] == 'Zero', "Last label should be the Zero class"
        self._names = single_cell.classes
        self._alpha = None
        self.config = config
        self.single_cell_data_missing = single_cell.isMissing

    @property
    def names(self) -> np.ndarray:
        """Returns the names of cell types."""
        assert self._names[-1] == 'Zero', "Last label should be the Zero class"
        return self._names

    @property
    def nK(self) -> int:
        """Returns the number of cell types."""
        return len(self.names)

    @property
    def alpha(self) -> np.ndarray:
        """Returns the alpha values for cell types."""
        return self._alpha

    @alpha.setter
    def alpha(self, val: np.ndarray):
        """Sets the alpha values for cell types."""
        self._alpha = val

    @property
    def pi_bar(self) -> np.ndarray:
        """Returns the pi bar values for cell types."""
        return self.alpha / self.alpha.sum()

    @property
    def logpi_bar(self) -> np.ndarray:
        """Returns the log pi bar values for cell types."""
        return scipy.special.psi(self.alpha) - scipy.special.psi(self.alpha.sum())

    @property
    def prior(self) -> np.ndarray:
        """Returns the prior probabilities for cell types."""
        return self.pi_bar

    @property
    def log_prior(self) -> np.ndarray:
        """Returns the log prior probabilities for cell types."""
        if self.single_cell_data_missing or self.config['cell_type_prior'] == 'weighted':
            return self.logpi_bar
        else:
            return np.log(self.prior)

    def size(self, cells: Cells) -> np.ndarray:
        """
        Calculates the size of each cell type, i.e., the number of cells in each type.

        Parameters:
            cells (Cells): Cells object containing cell data.

        Returns:
            np.array: Sizes of cell types.
        """
        return cells.classProb.sum(axis=0)

    def ini_prior(self):
        """Initializes the prior probabilities for cell types."""
        self.alpha = self.ini_alpha()

    def ini_alpha(self) -> np.ndarray:
        """
        Initializes the alpha values for cell types.

        Returns:
            np.array: Initialized alpha values.
        """
        ones = np.ones(self.nK - 1)
        return np.append(ones, sum(ones)).astype(np.float32)

    def calc_prior(self, total_counts, midpoint=60, steepness=0.15):
        """
        Calculate prior probabilities for all cell classes based on total gene counts.

        Uses a sigmoid function to assign higher Zero class probability to low-count cells
        and distributes remaining probability uniformly across non-Zero classes.

        Parameters:
        - total_counts: array of total counts per cell
        - midpoint: count value where Zero class probability = 0.5
        - steepness: how sharp the transition is (higher = sharper)

        Returns:
        - Array of shape (numCells, numClasses) with prior probabilities for all classes
        """
        zero_prob = 1 / (1 + np.exp((total_counts - midpoint) * steepness))

        # non-Zero cell types equally likely: uniform distribution across cell types (ex Zero class)
        cell_class_prob = (1-zero_prob) / (self.nK - 1)
        prob = np.column_stack([np.tile(cell_class_prob[:, None], self.nK-1), zero_prob])

        # assertion to verify probabilities sum to 1
        assert np.allclose(prob.sum(axis=1), 1.0), "Probabilities must sum to 1"

        return prob.astype(np.float32)

