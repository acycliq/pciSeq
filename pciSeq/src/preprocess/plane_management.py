"""
Plane management functionality for pciSeq.
Handles 3D data plane processing and spot filtering.
"""

from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import logging

plane_logger = logging.getLogger(__name__)


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
    mask_x = (spots.x >= 0) & (spots.x <= img_shape[2])
    mask_y = (spots.y >= 0) & (spots.y <= img_shape[1])
    mask_z = (spots.z_plane >= 0) & (spots.z_plane <= img_shape[0])
    return spots[mask_x & mask_y & mask_z]


def remove_planes(spots: pd.DataFrame,
                  coo: List[coo_matrix],
                  cfg: Dict) -> Tuple[pd.DataFrame, List[coo_matrix], int, pd.DataFrame]:
    """
    Remove specified planes from 3D data.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot data
    coo : List[coo_matrix]
        Label matrices
    cfg : Dict
        Configuration with exclude_planes

    Returns
    -------
    Tuple[pd.DataFrame, List[coo_matrix], int, pd.DataFrame]
        Processed spots, processed coo, minimum plane, removed cells
    """
    coo = label_image_remove_planes(coo, cfg)
    spots, min_plane = spots_remove_planes(spots, cfg)
    coo, removed = cells_remove_planes(coo, cfg)
    return spots, coo, min_plane, removed


def label_image_remove_planes(coo: List[coo_matrix], cfg: Dict) -> List[coo_matrix]:
    """Remove specified planes from label image."""
    arr = np.arange(len(coo))
    return [coo[d] for d in arr if d not in cfg['exclude_planes']]


def spots_remove_planes(spots: pd.DataFrame, cfg: Dict) -> Tuple[pd.DataFrame, int]:
    """
    Remove spots from excluded planes and adjust z coordinates.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot data
    cfg : Dict
        Configuration with exclude_planes

    Returns
    -------
    Tuple[pd.DataFrame, int]
        Processed spots and minimum plane number
    """
    int_z = np.floor(spots.z_plane)
    mask = [d not in cfg['exclude_planes'] for d in int_z]
    spots = spots[mask].copy()

    # Find first kept plane
    diff = np.diff(cfg['exclude_planes']) - 1
    if np.all(diff == 0):
        min_plane = max(cfg['exclude_planes']) + 1
    else:
        iLeft = list(diff > 0).index(True)
        min_plane = cfg['exclude_planes'][iLeft] + 1

    spots.loc[:, 'z_plane'] = spots.z_plane - min_plane
    return spots, min_plane


def cells_remove_planes(coo_list: List[coo_matrix],
                        cfg: Dict) -> Tuple[List[coo_matrix], pd.DataFrame]:
    """
    Remove cells that exist in only one frame.

    Parameters
    ----------
    coo_list : List[coo_matrix]
        Label matrices
    cfg : Dict
        Configuration

    Returns
    -------
    Tuple[List[coo_matrix], pd.DataFrame]
        Processed matrices and removed cell information
    """
    labels_per_frame = [np.unique(d.data) for d in coo_list]
    label_counts = np.bincount([d for labels in labels_per_frame for d in labels])
    single_page_labels = [d[0] for d in enumerate(label_counts) if d[1] == 1]

    removed_cells = []
    _frames = []

    for i, coo in enumerate(coo_list):
        s = set(coo.data).intersection(set(single_page_labels))
        for d in s:
            coo.data[coo.data == d] = 0
            coo.eliminate_zeros()
            removed_cells.append(d)
            _frames.append(i)

    if removed_cells:
        len_c = len(set(removed_cells))
        len_f = len(set(_frames))
        plane_logger.warning(
            f'Removed {len_c} cells that exist on single planes from {len_f} planes.'
        )

    removed_df = pd.DataFrame({
        'removed_cell_label': removed_cells,
        'frame_num': _frames,
        'comment': 'labels are the original labels from input segmentation masks'
    })

    return coo_list, removed_df
