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

import logging

my_logger = logging.getLogger(__name__)


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

    # make an extra column, the int of z_plane
    spots = spots.assign(plane_id=spots.z_plane.astype(np.int32))

    # make sure the index is consecutive
    assert spots.shape[0] == spots.index.max() + 1, "Spot indices are not consecutive"
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

    # Group by 'plane_id', apply the function, and reset the index so the result aligns with df.
    spots['label'] = (spots.groupby('plane_id')
                      .apply(inside_cell, coo)
                      .reset_index(level=0, drop=True).squeeze()
                      )
    return spots
