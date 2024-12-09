"""
Spot processing module for handling spot data transformations and assignments.
"""

from typing import List, Tuple
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from .label_processing import inside_cell
from ..core.utils.geometry import adjust_for_anisotropy
from .plane_management import remove_oob


def process_spots(spots: pd.DataFrame,
                  dimensions: Tuple[int, int, int],
                  voxel_size: Tuple[float, float, float]) -> pd.DataFrame:
    """
    Process spots by removing out-of-bounds and adjusting for anisotropy.

    Args:
        spots: DataFrame with spot coordinates
        dimensions: (n_planes, height, width) of image
        voxel_size: (x, y, z) voxel dimensions

    Returns:
        Processed spots DataFrame
    """
    spots = remove_oob(spots.copy(), dimensions)
    spots = adjust_for_anisotropy(spots, voxel_size)
    return spots


def assign_spot_labels(spots: pd.DataFrame, coo: List[coo_matrix]) -> pd.DataFrame:
    """
    Assign cell labels to spots based on their location.

    Args:
        spots: DataFrame with spot coordinates
        coo: List of sparse matrices containing cell labels

    Returns:
        Spots DataFrame with assigned labels
    """
    spots = spots.assign(label=np.zeros(spots.shape[0], dtype=np.uint32))

    for z in np.unique(spots.z_plane):
        spots_z = spots[spots.z_plane == z]
        inc = inside_cell(coo[int(z)].tocsr().astype(np.uint32), spots_z)
        spots.loc[spots.z_plane == z, 'label'] = inc

    return spots
