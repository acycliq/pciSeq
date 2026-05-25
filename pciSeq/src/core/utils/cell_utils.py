"""Utilities for processing and manipulating cell data."""
import sys
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any, Optional, Union
from numbers import Number
import math
from numba import njit, prange
import logging

# Configure logging
logger = logging.getLogger(__name__)


def read_image_objects(img_obj, cfg):
    meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
    relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

    # append 1 for the misreads
    relCellRadius = np.append(1, relCellRadius)

    InsideCellBonus = cfg['InsideCellBonus']
    if not InsideCellBonus:
        # This is more for clarity. The operation below will work fine even if InsideCellBonus is False
        InsideCellBonus = 0

    # if InsideCellBonus == 0 then CellAreaFactor will be equal to 1.0
    numer = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    denom = np.exp(-0.5) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    CellAreaFactor = numer / denom

    out = {
        'area_factor': CellAreaFactor.astype(np.float32), 'rel_radius': relCellRadius.astype(np.float32),
        'area': np.append(np.nan, img_obj.area.astype(np.uint32)),
        'x0': np.append(-sys.maxsize, img_obj.x0.values).astype(np.float32),
        'y0': np.append(-sys.maxsize, img_obj.y0.values).astype(np.float32),
        'z0': np.append(-sys.maxsize, img_obj.z0.values).astype(np.float32),
        'cell_label': np.append(0, img_obj.label.values).astype(np.uint32)
    }

    if 'old_label' in img_obj.columns:
        out['cell_label_old'] = np.append(0, img_obj.old_label.values).astype(np.uint32)
    # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
    # and will be used to get all the misreads. It was given the label=0 and some very small
    # negative coords

    return out, meanCellRadius.astype(np.float32)


def recover_original_labels(cellData: pd.DataFrame,
                            geneData: pd.DataFrame,
                            cellBoundaries: pd.DataFrame,
                            cellBoundaries_list: pd.DataFrame,
                            label_map: Optional[Dict[int, int]]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Restore original cell labels using label mapping.

    Args:
        cellData: Cell data DataFrame
        geneData: Gene data DataFrame
        label_map: Dictionary mapping new labels to original labels

    Returns:
        Tuple of (updated cellData, updated geneData)
    """
    if label_map is None:
        return cellData, geneData, cellBoundaries

    # Create reverse mapping
    reverse_map = {v: k for k, v in label_map.items()}

    # Update cell numbers
    cellData = cellData.assign(
        Cell_Num=cellData.Cell_Num.map(lambda x: reverse_map.get(x))
    )

    # Update gene data neighbors
    geneData = geneData.assign(
        neighbour=geneData.neighbour.map(lambda x: fetch_label(x, reverse_map)),
        neighbour_array=geneData.neighbour_array.map(lambda x: fetch_label(x, reverse_map))
    )

    cellBoundaries = cellBoundaries.assign(
        cell_id=cellBoundaries.cell_id.map(lambda x: reverse_map.get(x))
    )

    cellBoundaries_list = [d.assign(cell_id=d.cell_id.map(lambda x: reverse_map.get(x))) for d in cellBoundaries_list]

    logger.info("Restored original cell segmentation labels")
    return cellData, geneData, cellBoundaries, cellBoundaries_list


def fetch_label(x: Union[Number, List[Number]],
                d: Dict[int, int]) -> Union[int, List[int]]:
    """Fetch original label(s) from mapping dictionary.

    Args:
        x: Single label or list of labels
        d: Label mapping dictionary

    Returns:
        Original label(s)
    """
    x = [x] if isinstance(x, Number) else x
    out = [d[v] for v in x]
    return out[0] if len(out) == 1 else out


def keep_labels_unique(scdata: pd.DataFrame) -> pd.DataFrame:
    """Keep only highest count row for duplicate gene labels from the single cell reference data.

    Args:
        scdata: Single cell data DataFrame

    Returns:
        DataFrame with unique gene labels

    Notes:
        In single cell data you might find cases where two or more rows have
        the same gene label. In these cases keep the row with the highest
        total gene count.
    """
    # Add row total column
    scdata = scdata.assign(total=scdata.sum(axis=1))

    # Rank by gene label and total count, keep highest total
    scdata = (scdata.sort_values(['gene_name', 'total'],
                                 ascending=[True, False])
              .groupby('gene_name')
              .head(1))

    # Drop the total column and return
    return scdata.drop(['total'], axis=1)


@njit(parallel=True)
def create_circular_masks_large_scale(array_shape, centroids, radii):
    """
    Create circular masks on a 2D array for many centroids, with variable radii.
    """
    h, w = array_shape
    mask = np.zeros(array_shape, dtype=np.int8)

    # Process each centroid in parallel
    for i in prange(len(centroids)):
        center_x, center_y = centroids[i]
        radius = radii[i]
        radius_squared = radius ** 2

        # Calculate bounding box for this circle (with bounds checking)
        x_min = max(0, math.floor(center_x - radius))
        x_max = min(w, math.ceil(center_x + radius) + 1)
        y_min = max(0, math.floor(center_y - radius))
        y_max = min(h, math.ceil(center_y + radius) + 1)

        # Only process pixels within the bounding box
        for y in range(y_min, y_max):
            for x in range(x_min, x_max):
                dist_squared = (x - center_x) ** 2 + (y - center_y) ** 2
                if dist_squared <= radius_squared:
                    mask[y, x] = 1

    return mask


def create_circular_masks(array_shape, centroids, radii, chunk_size=1000):
    """
    Process centroids in chunks to reduce memory pressure, with flexible radius input.
    Parameters:
    - array_shape: Tuple (height, width) representing the shape of the array.
    - centroids: List of tuples [(x1, y1), (x2, y2), ...] representing the centers.
    - radii: Single number or list/array of radii corresponding to centroids.
    - chunk_size: Number of centroids to process in each batch.
    Returns:
    - A 2D array with the circular masks applied.
    """
    h, w = array_shape
    final_mask = np.zeros(array_shape, dtype=np.int8)

    # Convert inputs to numpy arrays
    centroids_array = np.array(centroids)

    # Handle different radii input types
    if np.isscalar(radii):
        # If a single number, create an array of repeated radii
        radii_array = np.full(len(centroids_array), radii)
    else:
        # Convert to numpy array
        radii_array = np.array(radii)

    # Validate input
    if len(centroids_array) != len(radii_array):
        raise ValueError("Number of centroids must match number of radii")

    # Process in chunks
    for i in range(0, len(centroids_array), chunk_size):
        chunk = centroids_array[i:i + chunk_size]
        chunk_radii = radii_array[i:i + chunk_size]
        chunk_mask = create_circular_masks_large_scale(array_shape, chunk, chunk_radii)
        final_mask = np.logical_or(final_mask, chunk_mask).astype(np.int8)

    return final_mask


def find_labels_by_plane_index(labels_dict, slice_idx):
    """
    Find all labels that appear in a specific slice.

    Args:
        labels_dict: Dictionary mapping label values to lists of slice indices
        slice_idx: The slice index to search for

    Returns:
        List of label values that appear in the specified slice
    """
    return [label for label, slices in labels_dict.items() if slice_idx in slices]


