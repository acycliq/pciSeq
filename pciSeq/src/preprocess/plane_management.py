"""
Plane management functionality for pciSeq.
Handles 3D data plane processing and spot filtering.
"""

from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import logging

logger = logging.getLogger(__name__)


def remove_oob(spots: pd.DataFrame, img_shape: List[int]) -> pd.DataFrame:
    """
    Remove out-of-bounds spots.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot coordinates
    img_shape : List[int]
        Image dimensions [z, y, x]

    Returns
    -------
    pd.DataFrame
        Filtered spots
    """
    mask_x = (spots.x >= 0) & (spots.x <= img_shape[2] - 1)
    mask_y = (spots.y >= 0) & (spots.y <= img_shape[1] - 1)
    mask_z = (spots.z_plane >= 0) & (spots.z_plane <= img_shape[0] - 1)
    return spots[mask_x & mask_y & mask_z]


def plane_quality_control(spots: pd.DataFrame,
                          coo: List[coo_matrix],
                          cfg: Dict) -> Tuple[pd.DataFrame, List[coo_matrix], pd.DataFrame]:
    """
    Perform quality control on 3D segmentation and spatial data.
    Removes single-plane cells.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot data
    coo : List[coo_matrix]
        Label matrices
    cfg : Dict
        Configuration

    Returns
    -------
    Tuple[pd.DataFrame, List[coo_matrix], pd.DataFrame]
        Processed spots, processed coo, removed cells
    """
    removed = pd.DataFrame()

    if cfg['remove_flat_cells']:
        coo, removed = remove_flat_cells(coo)
    return spots, coo, removed


def remove_flat_cells(coo_list: List[coo_matrix]) -> Tuple[List[coo_matrix], pd.DataFrame]:
    """
    Remove cells that exist in only one z-plane (segmentation artefacts).

    Edits the matrices in place and does the masking vectorised in-process. An
    earlier version deep-copied the whole stack and farmed the planes out to a
    multiprocessing Pool; both serialised the entire segmentation and dominated
    the runtime, while the actual work (flipping a handful of labels to zero) is
    tiny. We don't copy: nothing downstream needs the caller's coo_list pristine
    (process_labels also edits it in place right after), and on a big 3D stack a
    copy is wasted time and memory.

    Parameters
    ----------
    coo_list : List[coo_matrix]
        List of sparse matrices containing cell labels per z-plane.

    Returns
    -------
    Tuple[List[coo_matrix], pd.DataFrame]
        - The same matrices, edited in place, with single-plane cells removed.
        - DataFrame recording which cells were removed and from which planes.
    """
    # Fast path for empty input
    if not coo_list:
        return [], pd.DataFrame()

    # Validate input type
    if not all(isinstance(coo, coo_matrix) for coo in coo_list):
        raise ValueError("All elements in coo_list must be of type coo_matrix.")

    # 1: how many planes does each label appear in? Each plane's unique labels
    # concatenated, then a label whose total count is 1 lives in a single plane.
    per_plane_labels = np.concatenate([np.unique(coo.data) for coo in coo_list])
    labels, counts = np.unique(per_plane_labels, return_counts=True)
    single_page_labels = labels[counts == 1]

    # 2: zero those labels out, plane by plane, in place.
    removed_cells = []
    removed_planes = []
    if single_page_labels.size:
        for i, coo in enumerate(coo_list):
            mask = np.isin(coo.data, single_page_labels)
            if mask.any():
                removed_cells.extend(coo.data[mask].tolist())
                removed_planes.extend([i] * int(mask.sum()))
                coo.data[mask] = 0
                coo.eliminate_zeros()

    # 3: Log removal summary
    if removed_cells:
        logger.warning(
            f'Removed {len(set(removed_cells))} single-plane cells from {len(set(removed_planes))} planes.'
        )

    # 4: Create removal record
    removal_record = pd.DataFrame({
        'removed_cell_label': removed_cells,
        'frame_num': removed_planes,
        'comment': 'Original labels from segmentation masks'
    })

    return coo_list, removal_record
