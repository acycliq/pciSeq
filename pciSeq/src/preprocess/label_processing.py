"""
Label processing module for pciSeq preprocessing.
Handles all label-related operations including:
- Cell assignment
- Label normalization and mapping
- Label validation
- Label identification
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
import fastremap
import logging

label_processing_logger = logging.getLogger(__name__)


class CellLabelManager:
    """Manages cell label normalization and mapping."""

    def __init__(self, labels: np.ndarray):
        self.original_labels = labels
        self.label_map = None
        self._validate_labels()

    def _validate_labels(self) -> None:
        """Validate label integrity"""
        if not self._are_labels_sequential():
            label_processing_logger.warning('Non-sequential cell labels detected')

    def _are_labels_sequential(self) -> bool:
        """Check if labels are sequential"""
        unique_labels = np.unique(self.original_labels[self.original_labels > 0])
        return len(unique_labels) == unique_labels.max()

    def normalize_labels(self) -> Tuple[np.ndarray, Dict[int, int]]:
        """Create sequential labels and return mapping"""
        normalized_labels, label_map = fastremap.renumber(
            self.original_labels,
            in_place=False
        )
        self.label_map = label_map
        return normalized_labels, label_map

    def restore_original_labels(self, normalized_labels: np.ndarray) -> np.ndarray:
        """Restore original label numbering"""
        " This is not used, either find  way to use it or remove it."
        if self.label_map is None:
            return normalized_labels
        reverse_map = {v: k for k, v in self.label_map.items()}
        return fastremap.remap(normalized_labels, reverse_map, in_place=False)

    @classmethod
    def process_label_matrices(cls, coo_matrices: List[coo_matrix]) -> Tuple[List[coo_matrix], Optional[Dict]]:
        """
        Process cell label matrices to ensure sequential numbering.

        Args:
            coo_matrices: List of sparse matrices containing cell labels

        Returns:
            Tuple containing:
            - List of processed sparse matrices
            - Optional mapping dictionary if labels were renumbered
        """
        # Convert to dense array for processing
        arr_3d = np.stack([d.toarray() for d in coo_matrices])

        # Check and normalize labels if needed
        manager = cls(arr_3d)
        if not manager._are_labels_sequential():
            normalized_labels, label_map = manager.normalize_labels()
            label_processing_logger.warning('Labels have been renumbered for sequential labelling')
        else:
            normalized_labels, label_map = arr_3d, None

        # Convert back to sparse matrices
        processed_matrices = [coo_matrix(d) for d in normalized_labels]

        return processed_matrices, label_map


# def inside_cell(label_img: csr_matrix, spots: pd.DataFrame) -> np.ndarray:
#     """
#     Find which spots are inside cells.
#
#     Args:
#         label_img: Sparse matrix of cell labels
#         spots: DataFrame with spot coordinates
#
#     Returns:
#         Array of cell labels for each spot
#     """
#     x = spots.x.values
#     y = spots.y.values
#     label_dense = label_img.toarray()
#
#     # Get labels at spot coordinates
#     x = x.astype(int)
#     y = y.astype(int)
#     labels = label_dense[y, x]
#
#     return labels


def inside_cell(spots: pd.DataFrame, coo_list: List[coo_matrix]) -> pd.Series:
    """
    Compute labels for spots in a single plane group using the corresponding sparse matrix.

    Parameters
    ----------
    spots : pd.DataFrame
        DataFrame corresponding to a single plane group. Must have 'plane_id', 'x', and 'y' columns.
    coo_list : List[coo_matrix]
        List of sparse matrices containing cell labels.
    """
    unique_plane_ids = spots['plane_id'].unique()
    if len(unique_plane_ids) != 1:
        raise ValueError(f"Expected one unique plane_id per group, got: {unique_plane_ids}")
    plane_id = unique_plane_ids[0]

    # Convert the appropriate sparse matrix to CSR format.
    csr = coo_list[plane_id].tocsr()

    # Get the values at (y, x) positions and convert to a flattened 1D array.
    out = csr[spots['y'], spots['x']].A1

    # convert the list to a Series with the group's index. It needs to be a Series
    # or dataframe, so it will be properly aligned with the main spots dataframe
    return pd.Series(out, index=spots.index)


def get_unique_labels(coo_matrices: List[coo_matrix]) -> List[np.ndarray]:
    """
    Get unique labels from each image plane.

    Args:
        coo_matrices: List of sparse label matrices

    Returns:
        List of unique label arrays for each plane
    """
    return [np.unique(m.data) for m in coo_matrices if len(m.data) > 0]


def reorder_labels(labels: np.ndarray) -> Tuple[np.ndarray, Dict[int, int]]:
    """
    Reorder labels to be sequential.

    Args:
        labels: Array of cell labels

    Returns:
        Tuple containing:
        - Reordered labels array
        - Mapping from new to original labels
    """
    return fastremap.renumber(labels, in_place=False)
