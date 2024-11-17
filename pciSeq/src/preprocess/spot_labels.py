"""
Preprocessing module for pciSeq spatial transcriptomics data.

This module handles:
1. Cell label processing and validation
2. Spot-to-cell assignment
3. Cell boundary extraction
4. 3D data plane management
"""

from __future__ import annotations
from typing import List, Tuple, Dict, Optional, Union
import numpy as np
import pandas as pd
import skimage.measure as skmeas
import fastremap
from multiprocessing.dummy import Pool as ThreadPool
from scipy.sparse import coo_matrix, csr_matrix
from .cell_borders import extract_borders_dip, extract_borders
from ..core.utils import get_img_shape, adjust_for_anisotropy
import logging

spot_labels_logger = logging.getLogger(__name__)

# Type aliases for clarity
Labels = np.ndarray
CooList = List[coo_matrix]
SpotDF = pd.DataFrame
CellDF = pd.DataFrame


def inside_cell(label_image: Union[coo_matrix, csr_matrix, np.ndarray],
                spots: pd.DataFrame) -> np.ndarray:
    """
    Determine which spots lie within cell boundaries.

    Parameters
    ----------
    label_image : Union[coo_matrix, csr_matrix, np.ndarray]
        Image with cell labels
    spots : pd.DataFrame
        Spot coordinates

    Returns
    -------
    np.ndarray
        Cell labels for each spot
    """
    if isinstance(label_image, coo_matrix):
        label_image = label_image.tocsr()
    elif isinstance(label_image, np.ndarray):
        label_image = csr_matrix(label_image)
    elif isinstance(label_image, csr_matrix):
        pass
    else:
        raise TypeError('label_image should be of type "csr_matrix"')

    m = label_image[spots.y, spots.x]
    out = np.asarray(m, dtype=np.uint32)
    return out[0]


def get_unique_labels(masks: List[coo_matrix]) -> List[np.ndarray]:
    """
    Get unique cell labels from each mask.

    Parameters
    ----------
    masks : List[coo_matrix]
        List of label matrices

    Returns
    -------
    List[np.ndarray]
        Unique labels from each mask
    """
    pool = ThreadPool()
    out = pool.map(lambda mask: np.unique(mask.data), masks)
    pool.close()
    pool.join()
    return out


def reorder_labels(label_image: Union[np.ndarray, List[coo_matrix]]) -> Tuple[List[coo_matrix], Optional[pd.DataFrame]]:
    """
    Rearrange labels to form a sequence of integers.

    Parameters
    ----------
    label_image : Union[np.ndarray, List[coo_matrix]]
        Input label image

    Returns
    -------
    Tuple[List[coo_matrix], Optional[pd.DataFrame]]
        Reordered labels and mapping (if needed)
    """
    # Convert input to list of coo matrices if needed
    if isinstance(label_image, np.ndarray):
        if len(label_image.shape) not in {2, 3}:
            raise ValueError("Array must be 2D or 3D")
        label_image = [coo_matrix(d, dtype=np.uint32) for d in label_image]
    elif not (isinstance(label_image, list) and all(isinstance(d, coo_matrix) for d in label_image)):
        raise TypeError('Input must be array or list of coo matrices')

    # Get all unique labels
    labels = np.concatenate(get_unique_labels(label_image))

    # Check if reordering is needed
    if labels.max() != len(set(labels)):
        spot_labels_logger.info('Relabelling segmentation array to sequential integers...')
        labels = np.append(0, labels)  # include background
        _, idx = np.unique(labels, return_inverse=True)
        label_dict = dict(zip(labels, idx.astype(np.uint32)))

        if label_dict[0] != 0:
            raise ValueError("Background label (0) was not preserved")

        # Apply new labels
        out = []
        for plane in label_image:
            data = [label_dict[d] for d in plane.data]
            c = coo_matrix((data, (plane.row, plane.col)), shape=plane.shape)
            out.append(c)

        label_map = pd.DataFrame(label_dict.items(),
                                 columns=['old_label', 'new_label'],
                                 dtype=np.uint32)
        return out, label_map

    return label_image, None


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
        spot_labels_logger.warning(
            f'Removed {len_c} cells that exist on single planes from {len_f} planes.'
        )

    removed_df = pd.DataFrame({
        'removed_cell_label': removed_cells,
        'frame_num': _frames,
        'comment': 'labels are the original labels from input segmentation masks'
    })

    return coo_list, removed_df


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
