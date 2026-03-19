# Standard library imports
import logging
from typing import Dict

# Third party imports
import numpy as np
import scipy

from .singleCell import SingleCell
from .cells import Cells

logger = logging.getLogger(__name__)


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
        assert single_cell.classes[-1] == 'Zero', "Last cell class should be the Zero class"
        # Check that all classes except 'Zero' are in alphabetical order
        classes_without_zero = single_cell.classes[:-1]
        assert list(classes_without_zero) == sorted(classes_without_zero, key=str.lower), \
            "Cell type names (excluding 'Zero') must be in alphabetical order"
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
    def zero_weight(self) -> float:
        """Returns the fixed prior weight for the Zero class."""
        return self._initial_weights[-1]

    @property
    def pi_bar(self) -> np.ndarray:
        """
        Returns E[pi] for all classes, shape (nK,).

        Uniform mode: fixed initial weights.
        Weighted mode: Zero stays fixed, real classes from Dirichlet mean
            scaled to sum to (1 - zero_weight).
        """
        if self.config['cell_type_prior'] == 'uniform' and not self.single_cell_data_missing:
            return self._initial_weights

        # Dirichlet mean for real classes, scaled by (1 - zero_weight)
        alpha = self.alpha
        real_pi = (1 - self.zero_weight) * alpha / alpha.sum()
        return np.append(real_pi, self.zero_weight)

    @property
    def logpi_bar(self) -> np.ndarray:
        """
        Returns E[log pi] for all classes, shape (nK,).

        Uniform mode: log of fixed initial weights.
        Weighted mode: Zero gets log(zero_weight), real classes get
            log(1 - zero_weight) + psi(alpha_k) - psi(sum(alpha)).
        """
        if self.config['cell_type_prior'] == 'uniform' and not self.single_cell_data_missing:
            return np.log(self._initial_weights)

        # E[log Dir_k] for real classes, shifted by log(1 - zero_weight)
        alpha = self.alpha
        real_logpi = np.log(1 - self.zero_weight) + scipy.special.psi(alpha) - scipy.special.psi(alpha.sum())
        return np.append(real_logpi, np.log(self.zero_weight))

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
        self._initial_weights = self._compute_initial_weights()
        self.alpha = np.ones(self.nK - 1, dtype=np.float32)

    def _compute_initial_weights(self) -> np.ndarray:
        """
        Computes the initial probability vector for all classes from config.

        If cell_type_weights is None: flat 1/nK for all classes.
        If cell_type_weights is set: use the provided values, distribute remaining
        probability equally among unspecified classes.

        Returns:
            np.array: Probability vector of shape (nK,), sums to 1.
        """
        cfg_weights = self.config['cell_type_weights']
        if cfg_weights is None:
            return np.full(self.nK, 1.0 / self.nK, dtype=np.float64)

        specified = {}
        for key, val in cfg_weights.items():
            if key not in self.names:
                logger.warning(f"Cell type '{key}' in cell_type_weights not found in cell type names. Ignoring.")
            else:
                specified[key] = val

        specified_sum = sum(specified.values())
        if specified_sum > 1.0:
            logger.warning(f"cell_type_weights sum to {specified_sum} > 1.0, normalising.")
        remaining = max(0.0, 1.0 - specified_sum)
        unspecified = [name for name in self.names if name not in specified]
        default_weight = remaining / len(unspecified) if unspecified else 0.0

        weights = np.array([specified.get(name, default_weight) for name in self.names], dtype=np.float64)
        weights /= weights.sum()
        return weights
