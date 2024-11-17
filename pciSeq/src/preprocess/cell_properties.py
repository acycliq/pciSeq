"""
Cell property calculation functionality for pciSeq.
Handles measurement of cell characteristics and boundaries.
"""

from typing import List
import numpy as np
import pandas as pd
import skimage.measure as skmeas


def calculate_cell_properties(masks: np.ndarray, voxel_size: List[float]) -> pd.DataFrame:
    """
    Calculate cell properties from segmentation masks.

    Parameters
    ----------
    masks : np.ndarray
        3D array of cell labels
    voxel_size : List[float]
        Physical size of voxels [z, y, x]

    Returns
    -------
    pd.DataFrame
        Cell properties including position and size
    """
    scaling = [voxel_size[0] / voxel_size[0], voxel_size[1] / voxel_size[0], voxel_size[2] / voxel_size[0]]
    scaling = scaling[::-1]  # Convert to zyx order, same as the image

    properties = ['label', 'area', 'centroid', 'equivalent_diameter_area', 'bbox']
    props = skmeas.regionprops_table(
        label_image=masks,
        spacing=scaling,
        properties=properties
    )

    props_df = pd.DataFrame(props)
    props_df['mean_area_per_slice'] = (
            props_df['area'].values /
            (props_df['bbox-3'].values - props_df['bbox-0'].values)
    )

    props_df = props_df.rename(columns={
        "mean_area_per_slice": 'area',
        'area': 'volume',
        'centroid-0': 'z_cell',
        'centroid-1': 'y_cell',
        'centroid-2': 'x_cell'
    })

    props_df = props_df[['label', 'area', 'z_cell', 'y_cell', 'x_cell']]

    return props_df.astype({
        "label": np.uint32,
        "area": np.uint32,
        'z_cell': np.float32,
        'y_cell': np.float32,
        'x_cell': np.float32
    })
