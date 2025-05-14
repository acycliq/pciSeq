"""
Main preprocessing module for pciSeq spatial transcriptomics data.
Orchestrates the complete preprocessing pipeline.
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import logging
from .label_processing import CellLabelManager, get_unique_labels
from .spot_processing import process_spots, assign_spot_labels
from .utils import log_data_summary
from .plane_management import plane_quality_control
from .cell_processing import calculate_cell_properties
from ..core.utils.geometry import get_img_shape
from ..core.utils.cell_utils import find_labels
from .cell_processing import extract_borders

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
    label_map : Optional[Dict]
        Label remapping if labels were reordered
    """
    # Perform quality control on 3D data
    if cfg['is3D']:
        spots, coo, min_plane, removed = plane_quality_control(spots, coo, cfg)
        if removed.shape[0] > 0:
            removed.frame_num = removed.frame_num + min_plane

    # Process label matrices
    coo, label_map = CellLabelManager.process_label_matrices(coo)
    cfg['label_map'] = label_map  # Dont quite like it here, need to make it more transparent!!

    img_dim = {'n_planes': len(coo),
               'w': coo[0].shape[1],
               'h': coo[0].shape[0]}
    cfg['img_dim'] = img_dim

    # Process spots
    dimensions = get_img_shape(coo)
    spots = process_spots(spots, dimensions, cfg['voxel_size'])
    log_data_summary(spots, coo, dimensions)
    spots = assign_spot_labels(spots, coo)

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
    spot_labels_logger.info("find_labels start")
    labels_dict = find_labels(masks)
    spot_labels_logger.info("find_labels end")
    labels_df = pd.DataFrame(labels_dict.items(), columns=['index', 'values']).set_index('index')

    cells = props_df.rename(columns={'x_cell': 'x0', 'y_cell': 'y0', 'z_cell': 'z0'})

    assert np.all(labels_df.index.values == cells.label.values)
    cells = cells.merge(labels_df, how='left', left_on='label', right_on='index')
    processed_spots = spots[['x', 'y', 'z', 'plane_id', 'label', 'gene_name', 'score', 'intensity']].rename_axis('spot_id')

    return cells, cell_boundaries, processed_spots, label_map
