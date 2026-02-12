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
        assert single_cell.classes[-1] == 'Zero', "Last cell class should be the Zero class"
        # Check that all classes except 'Zero' are in alphabetical order
        classes_without_zero = single_cell.classes[:-1]
        assert list(classes_without_zero) == sorted(classes_without_zero, key=str.lower), \
            "Cell type names (excluding 'Zero') must be in alphabetical order"
        self._names = single_cell.classes
        self._alpha = None
        self._prior = None
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
        return self._prior

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
        cfg_weights = self.config['cell_type_weights']
        if cfg_weights:
            # Extract default value if provided, otherwise use 1
            default_weight = cfg_weights.get('default', 1)

            # Initialize all cell types with the default weight
            weights_dict = {name: default_weight for name in self.names}

            for key in cfg_weights:
                if key == 'default':
                    # Skip the 'default' key as it's not a cell type
                    continue
                elif key not in self.names:
                    cellType_logger.warning(f"Cell type '{key}' in cell_type_weights not found in cell type names. Ignoring.")
                else:
                    # Override with provided weights where applicable
                    weights_dict[key] = cfg_weights[key]

            # Preserve the 'Zero = sum(others)' behavior unless explicitly overridden
            if 'Zero' not in cfg_weights:
                vals = list(weights_dict.values())
                weights_dict["Zero"] = np.sum(vals[:-1])

            # get the values from the dict as a numpy array
            out = np.array(list(weights_dict.values()), dtype=np.float32)
        else:
            ones = np.ones(self.nK - 1)
            out = np.append(ones, sum(ones)).astype(np.float32)

        return out


    def ini_prior_v2(self, cell_centroids):
        weight_dict = {
            "CA1": {"016 CA1-ProS Glut": 0.6},
            "CA2": {"025 CA2-FC-IG Glut": 0.6},
            "CA3": {"017 CA3 Glut": 0.6},
            "DG":  {"037 DG Glut": 0.3, "038 DG-PIR Ex IMN": 0.3},
        }
        region_labels = self.mask_cells(cell_centroids)
        classes = self.names
        class_to_idx = {c: j for j, c in enumerate(classes)}
        out = np.zeros((len(region_labels), self.nK))

        for i, region in enumerate(region_labels):
            class_weight = weight_dict.get(region,  {"Zero": 0.5})
            total = sum(class_weight.values())
            remaining = 1.0 - total

            # uniform fill across non-fixed classes
            fixed_indices = {class_to_idx[c] for c in class_weight}
            other_count = self.nK - len(fixed_indices)
            if other_count > 0:
                out[i, :] = remaining / other_count

            # overwrite fixed classes
            for cls_name, prob in class_weight.items():
                out[i, class_to_idx[cls_name]] = prob
        self._prior = out

    def mask_cells(self, centroids):
        from shapely.geometry import Polygon
        import shapely
        import pandas as pd

        centroid_points = shapely.points(centroids[['x', 'y']].values)

        # Load bounding box polygons
        ca1_bbox = pd.read_csv('./silver_metadata/region_boundaries/ca1_bbox.csv')
        ca2_bbox = pd.read_csv('./silver_metadata/region_boundaries/ca2_bbox.csv')
        ca3_bbox = pd.read_csv('./silver_metadata/region_boundaries/ca3_bbox.csv')
        dg_bbox = pd.read_csv('./silver_metadata/region_boundaries/dg_bbox.csv')

        ca1_polygon = Polygon(ca1_bbox.values)
        ca2_polygon = Polygon(ca2_bbox.values)
        ca3_polygon = Polygon(ca3_bbox.values)
        dg_polygon = Polygon(dg_bbox.values)

        # Add spatial containment columns
        in_ca1=shapely.contains(ca1_polygon, centroid_points)
        in_ca2=shapely.contains(ca2_polygon, centroid_points)
        in_ca3=shapely.contains(ca3_polygon, centroid_points)
        in_dg=shapely.contains(dg_polygon, centroid_points)

        out = np.full(len(centroid_points), "Other", dtype=object)
        out[in_ca1] = "CA1"
        out[in_ca2] = "CA2"
        out[in_ca3] = "CA3"
        out[in_dg] = "DG"

        return out
