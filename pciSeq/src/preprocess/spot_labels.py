"""
Main preprocessing module for pciSeq spatial transcriptomics data.
Orchestrates the complete preprocessing pipeline.
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import fastremap
import logging
from .label_processing import inside_cell, get_unique_labels, reorder_labels
from .plane_management import remove_oob, remove_planes
from .cell_properties import calculate_cell_properties
from ..core.utils import get_img_shape, adjust_for_anisotropy
from .cell_borders import extract_borders

spot_labels_logger = logging.getLogger(__name__)


def stage_data(spots: pd.DataFrame,
               coo: List[coo_matrix],
               cfg: Dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Optional[Dict]]:
    """
    Process spots and label images for cell typing analysis.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot data with columns: ['gene_name', 'x', 'y', 'z_plane']
    coo : List[coo_matrix]
        List of sparse matrices containing cell segmentation
    cfg : Dict
        Configuration dictionary with processing parameters

    Returns
    -------
    cells : pd.DataFrame
        Cell properties including position and size
    cell_boundaries : pd.DataFrame
        Cell boundary coordinates
    processed_spots : pd.DataFrame
        Processed spots with cell assignments
    remapping : Optional[Dict]
        Label remapping if labels were reordered
    """
    remapping = None

    # Handle 3D data plane exclusion
    if cfg['is3D'] and cfg['exclude_planes'] is not None:
        spots, coo, min_plane, removed = remove_planes(spots, coo, cfg)
        if removed.shape[0] > 0:
            removed.frame_num = removed.frame_num + min_plane

    # Check for skipped labels
    arr_3d = np.stack([d.toarray() for d in coo])
    if arr_3d.max() != sum(np.unique(arr_3d) > 0):
        spot_labels_logger.warning('Skipped labels detected in the label image.')
        arr_3d, remapping = fastremap.renumber(arr_3d, in_place=True)
        coo = [coo_matrix(d) for d in arr_3d]
        del arr_3d

    # Process spots
    [n, h, w] = get_img_shape(coo)
    spots = remove_oob(spots.copy(), [n, h, w])
    spots = adjust_for_anisotropy(spots, cfg['voxel_size'])

    # Log data summary
    spot_labels_logger.info(f'Number of spots passed-in: {spots.shape[0]}')
    spot_labels_logger.info(f'Number of segmented cells: {max([d.data.max() for d in coo if len(d.data) > 0])}')
    if n == 1:
        spot_labels_logger.info(f'Image dimensions: {w}px × {h}px')
    else:
        spot_labels_logger.info(f'Image dimensions: {n} planes, {w}px × {h}px')

    # Assign spots to cells
    spots = spots.assign(label=np.zeros(spots.shape[0], dtype=np.uint32))
    for z in np.unique(spots.z_plane):
        spots_z = spots[spots.z_plane == z]
        inc = inside_cell(coo[int(z)].tocsr().astype(np.uint32), spots_z)
        spots.loc[spots.z_plane == z, 'label'] = inc

    # Calculate cell properties
    masks = np.stack([d.toarray().astype(np.uint32) for d in coo])
    props_df = calculate_cell_properties(masks, cfg['voxel_size'])

    # Get cell boundaries
    mid_plane = int(np.floor(len(coo) / 2))
    cell_boundaries = extract_borders(coo[mid_plane].toarray().astype(np.uint32))
    cell_boundaries = cell_boundaries.rename(columns={'label': 'cell_id'})

    # Validate results
    assert props_df.shape[0] == len(set(np.concatenate(get_unique_labels(coo))))
    assert set(spots.label[spots.label > 0]) <= set(props_df.label)

    # Prepare output
    cells = props_df.rename(columns={'x_cell': 'x0', 'y_cell': 'y0', 'z_cell': 'z0'})
    processed_spots = spots[['x', 'y', 'z', 'label', 'gene_name']].rename_axis('spot_id')

    return cells, cell_boundaries, processed_spots, remapping
