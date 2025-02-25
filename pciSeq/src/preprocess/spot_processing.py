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

    my_logger.info('my_inside starts')

    # Group by 'plane_id', apply the function, and reset the index so the result aligns with df.
    spots['label'] = (spots.groupby('plane_id')
                      .apply(inside_cell, coo)
                      .reset_index(level=0, drop=True)
                      )
    my_logger.info('my_inside finished')

    # my_logger.info('inside_cell loop starts')
    # for z in np.unique(spots.z_plane):
    #     spots_z = spots[spots.z_plane == z]
    #     inc = inside_cell(coo[int(z)].tocsr().astype(np.uint32), spots_z)
    #     spots.loc[spots.z_plane == z, 'label'] = inc
    #
    # my_logger.info('inside_cell loop finished')
    return spots

# def my_inside(spots, coo):
#     pid = set(spots.plane_id)
#     print(pid)
#     assert len(pid) == 1
#     pid = pid.pop()
#     csr = coo[pid].tocsr()
#     out = csr[spots.y, spots.x]
#     out = out.tolist()[0]
#
#     # convert the list to a Series with the group's index
#     return pd.Series(out, index=spots.index)


# def my_inside(spots: pd.DataFrame, coo_list: List[coo_matrix]) -> pd.Series:
#     """
#     Compute labels for spots in a single plane group using the corresponding sparse matrix.
#
#     Parameters
#     ----------
#     spots : pd.DataFrame
#         DataFrame corresponding to a single plane group. Must have 'plane_id', 'x', and 'y' columns.
#     coo_list : List[coo_matrix]
#         List of sparse matrices containing cell labels.
#     """
#     unique_plane_ids = spots['plane_id'].unique()
#     if len(unique_plane_ids) != 1:
#         raise ValueError(f"Expected one unique plane_id per group, got: {unique_plane_ids}")
#     plane_id = unique_plane_ids[0]
#
#     # Convert the appropriate sparse matrix to CSR format.
#     csr = coo_list[plane_id].tocsr()
#
#     # Get the values at (y, x) positions and convert to a flattened 1D array.
#     out = csr[spots['y'], spots['x']].A1
#
#     # convert the list to a Series with the group's index. It needs to be a Series
#     # or dataframe, so it will be properly aligned with the main spots dataframe
#     return pd.Series(out, index=spots.index)
