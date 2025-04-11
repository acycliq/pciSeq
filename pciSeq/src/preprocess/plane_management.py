"""
Plane management functionality for pciSeq.
Handles 3D data plane processing and spot filtering.
"""

from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from collections import defaultdict
from copy import deepcopy
from multiprocessing import Pool, cpu_count
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
    mask_x = (spots.x >= 0) & (spots.x <= img_shape[2] - 1)
    mask_y = (spots.y >= 0) & (spots.y <= img_shape[1] - 1)
    mask_z = (spots.z_plane >= 0) & (spots.z_plane <= img_shape[0] - 1)
    return spots[mask_x & mask_y & mask_z]


def plane_quality_control(spots: pd.DataFrame,
                          coo: List[coo_matrix],
                          cfg: Dict) -> Tuple[pd.DataFrame, List[coo_matrix], int, pd.DataFrame]:
    """
    Perform quality control on 3D segmentation and spatial data.
    Handles plane exclusion and removes single-plane cells.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot data
    coo : List[coo_matrix]
        Label matrices
    cfg : Dict
        Configuration with optional exclude_planes

    Returns
    -------
    Tuple[pd.DataFrame, List[coo_matrix], int, pd.DataFrame]
        Processed spots, processed coo, minimum plane, removed cells
    """
    min_plane = 0
    removed = pd.DataFrame()
    if cfg['exclude_planes'] is not None:
        coo = label_image_remove_planes(coo, cfg)
        spots, min_plane = spots_remove_planes(spots, cfg)

    if cfg['remove_flat_cells']:
        coo, removed = remove_flat_cells_par(coo)
    return spots, coo, min_plane, removed


def label_image_remove_planes(coo: List[coo_matrix], cfg: Dict) -> List[coo_matrix]:
    """Remove specified planes from label image."""
    arr = np.arange(len(coo))
    return [coo[d] for d in arr if d not in cfg['exclude_planes']]


def spots_remove_planes(spots: pd.DataFrame, cfg: Dict) -> Tuple[pd.DataFrame, int]:
    """
    Remove spots from excluded planes and adjust z coordinates.
    !!!!!!! MUST BE REVIEWED !!!!!

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
    if cfg['exclude_planes']:
        diff = np.diff(cfg['exclude_planes']) - 1
        if np.all(diff == 0):
            min_plane = max(cfg['exclude_planes']) + 1
        else:
            iLeft = list(diff > 0).index(True)
            min_plane = cfg['exclude_planes'][iLeft] + 1

        spots.loc[:, 'z_plane'] = spots.z_plane - min_plane
    else:
        min_plane = 0

    return spots, min_plane


def remove_flat_cells(coo_list: List[coo_matrix]) -> Tuple[List[coo_matrix], pd.DataFrame]:
    """
    Remove cells that exist in only one plane

    Parameters
    ----------
    coo_list : List[coo_matrix]
        List of sparse matrices containing cell labels per z-plane

    Returns
    -------
    Tuple[List[coo_matrix], pd.DataFrame]
        - Modified matrices with single-plane cells removed
        - DataFrame recording which cells were removed and from which planes
    """
    # Fast path for empty input
    if not coo_list:
        return [], pd.DataFrame()

    # 1: Identify single-plane cells
    # 1.1: Get all unique labels present in each plane
    labels_per_frame = [np.unique(d.data) for d in coo_list]
    # 1.2: Count how many times each label appears across all planes
    label_counts = np.bincount([d for labels in labels_per_frame for d in labels])
    # 1.3: Get the labels that appear in only one plane
    single_page_labels = set(d for d, count in enumerate(label_counts) if count == 1)

    # 2: Process each plane and track removals
    removed_cells = []
    removed_planes = []

    for i, coo in enumerate(coo_list):
        # Find intersection of current plane's labels with single-plane labels
        intersected_labels = set(coo.data).intersection(single_page_labels)
        for label in intersected_labels:
            # set all occurrences of the current label to zero
            coo.data[coo.data == label] = 0
            coo.eliminate_zeros()
            # record keeping
            removed_cells.append(label)
            removed_planes.append(i)

    # 3: Log removal summary
    if removed_cells:
        plane_logger.warning(
            f'Removed {len(set(removed_cells))} single-plane cells from {len(set(removed_planes))} planes.'
        )

    # Step 4: Create removal record
    removal_record = pd.DataFrame({
        'removed_cell_label': removed_cells,
        'frame_num': removed_planes,
        'comment': 'Original labels from segmentation masks'
    })

    return coo_list, removal_record


def process_plane(args):
    """
    Helper function to process a single plane in parallel.

    Parameters
    ----------
    args : Tuple[int, coo_matrix, set]
        A tuple containing:
        - Index of the plane (int)
        - The sparse matrix (coo_matrix)
        - Set of single-plane labels to remove (set)

    Returns
    -------
    Tuple[int, coo_matrix, List[int]]
        - Index of the plane (int)
        - Modified sparse matrix (coo_matrix)
        - List of removed cell labels (List[int])
    """
    i, coo, single_page_labels = args
    # Find intersection of current plane's labels with single-plane labels
    mask = np.isin(coo.data, list(single_page_labels))
    removed_cells = coo.data[mask].tolist() if np.any(mask) else []
    # Remove single-plane cells
    coo.data[mask] = 0
    coo.eliminate_zeros()
    return i, coo, removed_cells


def remove_flat_cells_par(coo_list: List[coo_matrix]) -> Tuple[List[coo_matrix], pd.DataFrame]:
    """
    Remove cells that exist in only one plane. Parallelised version of remove_flat_cells

    Parameters
    ----------
    coo_list : List[coo_matrix]
        List of sparse matrices containing cell labels per z-plane.

    Returns
    -------
    Tuple[List[coo_matrix], pd.DataFrame]
        - Modified matrices with single-plane cells removed.
        - DataFrame recording which cells were removed and from which planes.
    """
    # Fast path for empty input
    if not coo_list:
        return [], pd.DataFrame()

    # Validate input type
    if not all(isinstance(coo, coo_matrix) for coo in coo_list):
        raise ValueError("All elements in coo_list must be of type coo_matrix.")

    # Make a deep copy of the input to avoid in-place modification
    coo_list = deepcopy(coo_list)

    # 1: Identify single-plane cells
    # Use a dictionary to count occurrences of each label across all planes
    label_counts = defaultdict(int)
    for coo in coo_list:
        unique_labels = np.unique(coo.data)
        for label in unique_labels:
            label_counts[label] += 1

    # Get the labels that appear in only one plane
    single_page_labels = {label for label, count in label_counts.items() if count == 1}

    # 2: Process each plane in parallel
    removed_cells = []
    removed_planes = []

    # Prepare arguments for parallel processing
    args = [(i, coo, single_page_labels) for i, coo in enumerate(coo_list)]

    # Use multiprocessing to process planes in parallel
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_plane, args)

    # Reconstruct the modified coo_list and track removals
    for i, coo, cells in results:
        coo_list[i] = coo
        if cells:
            removed_cells.extend(cells)
            removed_planes.extend([i] * len(cells))

    # 3: Log removal summary
    if removed_cells:
        plane_logger.warning(
            f'Removed {len(set(removed_cells))} single-plane cells from {len(set(removed_planes))} planes.'
        )

    # 4: Create removal record
    removal_record = pd.DataFrame({
        'removed_cell_label': removed_cells,
        'frame_num': removed_planes,
        'comment': 'Original labels from segmentation masks'
    })

    return coo_list, removal_record
