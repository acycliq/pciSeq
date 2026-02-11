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


def process_labels(coo_list: List[coo_matrix]) -> Tuple[List[coo_matrix], Optional[Dict]]:
    """
    Ensure cell labels are sequential (1..N) across all planes.
    Works directly on sparse .data arrays — no dense conversion needed.

    Args:
        coo_list: List of sparse matrices containing cell labels

    Returns:
        Tuple containing:
        - List of processed sparse matrices (modified in-place)
        - Optional mapping dictionary if labels were renumbered
    """
    all_data = np.concatenate([coo.data for coo in coo_list if coo.nnz > 0])

    unique_labels = fastremap.unique(all_data)
    unique_labels = unique_labels[unique_labels > 0]
    is_sequential = len(unique_labels) == unique_labels.max()

    if not is_sequential:
        label_processing_logger.warning('Non-sequential cell labels detected')
        _, label_map = fastremap.renumber(
            all_data, in_place=False, preserve_zero=True
        )

        # Remap each sparse matrix's .data array in-place. Since numpy arrays
        # are mutable, this modifies the caller's data without copying.
        for coo in coo_list:
            if coo.nnz > 0:
                fastremap.remap(coo.data, label_map, in_place=True)
        label_processing_logger.warning('Labels have been renumbered for sequential labelling')
    else:
        label_map = None

    return coo_list, label_map


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


