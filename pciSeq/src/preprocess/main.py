"""
Main preprocessing module for pciSeq spatial transcriptomics data.
Orchestrates the complete preprocessing pipeline.
"""

from typing import List, Tuple, Dict
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import logging
from joblib import Parallel, delayed
from .label_processing import process_labels
from .spot_processing import process_spots, assign_spot_labels
from .utils import log_data_summary
from .plane_management import plane_quality_control
from .cell_processing import calculate_cell_properties
from ..core.utils.geometry import get_img_shape
from .cell_processing import extract_borders

logger = logging.getLogger(__name__)


def _process_plane_borders(i: int, coo_plane: coo_matrix) -> pd.DataFrame:
    """
    Helper function to extract borders for a single plane.

    Parameters
    ----------
    i : int
        Plane index
    coo_plane : coo_matrix
        Sparse matrix for this plane

    Returns
    -------
    pd.DataFrame
        Borders dataframe with plane_id column
    """
    temp = extract_borders(coo_plane.toarray().astype(np.uint32))
    temp = temp.rename(columns={'label': 'cell_id'})
    temp.insert(0, 'plane_id', i)
    return temp


def _extract_all_borders(coo: List[coo_matrix]) -> Tuple[pd.DataFrame, List[pd.DataFrame]]:
    """Extract borders for all planes. Runs in a background thread."""
    logger.info("Border extraction started...")
    mid_plane = len(coo) // 2
    cell_boundaries_list = Parallel(n_jobs=-1, backend='loky')(
        delayed(_process_plane_borders)(i, d) for i, d in enumerate(coo)
    )
    cell_boundaries = cell_boundaries_list[mid_plane]
    logger.info("Border extraction complete")
    return cell_boundaries, cell_boundaries_list


def stage_data(spots: pd.DataFrame,
               coo: List[coo_matrix],
               cfg: Dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
    borders_future : Future
        Future resolving to (cell_boundaries, cell_boundaries_list).
        Border extraction runs in the background and only blocks when .result() is called.
    processed_spots : pd.DataFrame
        Processed spots with cell assignments

    Note
    ----
    The label remapping (label_map) and image dimensions (img_dim) are written
    into cfg as runtime state, the same way Config.set_runtime_attrs adds is3D.
    """
    # Perform quality control on 3D data
    if cfg['is3D']:
        # the removal record is logged inside plane_quality_control; not needed here
        spots, coo, _ = plane_quality_control(spots, coo, cfg)

    # Process label matrices
    coo, label_map = process_labels(coo)
    # runtime-derived run state, kept in cfg alongside is3D (see Config.set_runtime_attrs)
    cfg['label_map'] = label_map

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
    props_df = calculate_cell_properties(coo, cfg['voxel_size'])

    # Launch border extraction in the background, not needed until results are saved
    logger.info(f"Submitting border extraction for {len(coo)} planes (background task)...")
    executor = ThreadPoolExecutor(max_workers=1)
    borders_future = executor.submit(_extract_all_borders, coo)
    executor.shutdown(wait=False)

    # Validate results. Explicit raises (not asserts) so the checks still run
    # under python -O; on real data we always want these integrity checks active.
    # Both checks use cheap numpy ops, not python set() over millions of rows.
    #
    # n_labels = how many distinct cell labels are in the stack. We mark them in a
    # presence array (a scatter, O(nnz)) instead of np.unique per plane (a sort per
    # plane), which is a lot cheaper on a big 3D stack.
    max_label = max((int(m.data.max()) for m in coo if m.nnz), default=0)
    present = np.zeros(max_label + 1, dtype=bool)
    for m in coo:
        if m.nnz:
            present[m.data] = True
    n_labels = int(present[1:].sum())  # exclude background (label 0)
    if props_df.shape[0] != n_labels:
        raise RuntimeError(
            f"cell property rows ({props_df.shape[0]}) do not match the number of cell labels ({n_labels})"
        )

    # every spot must be assigned to a cell label that has computed properties
    spot_labels = np.unique(spots.label.values)
    spot_labels = spot_labels[spot_labels > 0]
    if not np.isin(spot_labels, props_df.label.values).all():
        raise RuntimeError("some spots are assigned to cell labels with no computed properties")

    cells = props_df.rename(columns={'x_cell': 'x0', 'y_cell': 'y0', 'z_cell': 'z0'})
    processed_spots = spots[['x', 'y', 'z', 'plane_id', 'label', 'gene_name', 'score', 'intensity']].rename_axis('spot_id')

    return cells, borders_future, processed_spots
