"""
Cell and Spot Label Processing Module for pciSeq

This module provides core functionality for processing and analyzing spatial transcriptomics data,
specifically handling the relationship between cell segmentation and RNA spot detection.

Key Functions:
-------------
- inside_cell: Maps RNA spots to their containing cells
- reorder_labels: Normalizes cell labels to sequential integers
- stage_data: Main processing pipeline that integrates all functionality

Data Processing Steps:
--------------------
1. Cell Label Processing:
   - Validates and normalizes cell segmentation labels
   - Ensures sequential integer labeling
   - Computes cell properties (centroids, areas)

2. Spot Assignment:
   - Maps each RNA spot to its containing cell
   - Validates spot coordinates against image boundaries
   - Links spots with cell metadata

3. Boundary Processing:
   - Extracts and validates cell boundaries
   - Ensures consistency between properties and boundaries

"""

from dataclasses import dataclass


import numpy as np
import pandas as pd
from typing import Tuple, Union
import skimage.measure as skmeas
from scipy.sparse import coo_matrix, csr_matrix, spmatrix
from pciSeq.src.preprocess.cell_borders import extract_borders
import logging


spot_labels_logger = logging.getLogger(__name__)


def inside_cell(label_image: Union[spmatrix, np.ndarray],
                spots: pd.DataFrame) -> np.ndarray:
    """
    Determine which cell contains each spot.

    Args:
        label_image: Cell segmentation mask (sparse matrix or array)
        spots: DataFrame with spot coordinates ('x' and 'y' columns)

    Returns:
        Array of cell labels for each spot

    Raises:
        TypeError: If label_image is not of supported type
    """
    if isinstance(label_image, np.ndarray):
        label_image = csr_matrix(label_image)
    elif isinstance(label_image, coo_matrix):
        label_image = label_image.tocsr()
    elif not isinstance(label_image, csr_matrix):
        raise TypeError('label_image must be ndarray, coo_matrix, or csr_matrix')

    return np.asarray(label_image[spots.y, spots.x], dtype=np.uint32)[0]


def remap_labels(coo: coo_matrix) -> coo_matrix:
    """
    Randomly reshuffle label assignments (for testing/debugging).

    Args:
        coo: Input label matrix

    Returns:
        Matrix with randomly remapped labels
    """
    coo_max = coo.data.max()
    original_labels = 1 + np.arange(coo_max)
    new_labels = original_labels.copy()
    np.random.shuffle(new_labels)

    label_map = dict(zip(original_labels, new_labels))
    new_data = np.array([label_map[x] for x in coo.data]).astype(np.uint64)

    return coo_matrix((new_data, (coo.row, coo.col)), shape=coo.shape)


def reorder_labels(coo: coo_matrix) -> Tuple[coo_matrix, pd.DataFrame]:
    """
    Normalize labels to be sequential integers starting from 1.

    Args:
        coo: Sparse matrix containing cell labels

    Returns:
        Tuple containing:
            - Relabeled sparse matrix
            - DataFrame mapping old labels to new labels
    """
    label_image = coo.toarray()
    unique_labels, idx = np.unique(label_image.flatten(), return_inverse=True)

    label_map = pd.DataFrame({
        'old_label': unique_labels,
        'new_label': np.arange(len(unique_labels))
    }, dtype=np.uint32).sort_values(by='old_label', ignore_index=True)

    return coo_matrix(idx.reshape(label_image.shape)), label_map


def stage_data(spots: pd.DataFrame, coo: coo_matrix) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Process spot and cell segmentation data to establish spot-cell relationships.

    Args:
        spots: DataFrame with columns ['x', 'y', 'Gene']
        coo: Sparse matrix containing cell segmentation labels

    Returns:
        Tuple containing:
            - cells: DataFrame with cell properties
            - boundaries: DataFrame with cell boundary coordinates
            - spots: DataFrame with spot locations and cell assignments

    Raises:
        ValueError: If required columns are missing or data validation fails
    """
    # Validate inputs
    required_columns = {'x', 'y', 'Gene'}
    missing_columns = required_columns - set(spots.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns in spots DataFrame: {missing_columns}")

    spot_labels_logger.info(f'Number of spots passed-in: {spots.shape[0]}')
    spot_labels_logger.info(f'Number of segmented cells: {len(set(coo.data))}')
    spot_labels_logger.info(
        f'Segmentation array implies that image has width: {coo.shape[1]}px and height: {coo.shape[0]}px')

    # Normalize labels if needed
    label_map = None
    unique_labels = set(coo.data)
    max_label = coo.data.max()

    if coo.data.max() != len(set(coo.data)):
        spot_labels_logger.info(
            f'Detected non-sequential cell labels: found {len(unique_labels)} unique labels '
            f'with maximum value of {max_label}. Normalizing labels to range [1, {len(unique_labels)}]'
        )
        coo, label_map = reorder_labels(coo)

    # Filter spots to image bounds
    mask_x = (spots.x >= 0) & (spots.x <= coo.shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= coo.shape[0])
    spots = spots[mask_x & mask_y].copy()

    # Assign spots to cells
    spots = spots.assign(label=inside_cell(coo, spots))

    # Calculate cell properties
    props = skmeas.regionprops_table(
        coo.toarray().astype(np.int32),
        properties=['label', 'area', 'centroid']
    )
    props_df = pd.DataFrame(props).rename(
        columns={'centroid-0': 'y_cell', 'centroid-1': 'x_cell'}
    )

    # Apply label mapping if exists
    if label_map is not None:
        props_df = pd.merge(
            props_df,
            label_map,
            left_on='label',
            right_on='new_label',
            how='left'
        ).drop(['new_label'], axis=1)

    # Set datatypes
    props_df = props_df.astype({
        "label": np.uint32,
        "area": np.uint32,
        'y_cell': np.float32,
        'x_cell': np.float32
    })

    # Extract cell boundaries
    cell_boundaries = extract_borders(coo.toarray().astype(np.uint32))

    # Validate cell data consistency
    if not (props_df.shape[0] == cell_boundaries.shape[0] == np.unique(coo.data).shape[0]):
        raise ValueError("Inconsistency detected between cell properties and boundaries")

    # Ensure spots are assigned to valid cells
    if not set(spots.label[spots.label > 0]).issubset(set(props_df.label)):
        raise ValueError("Spots assigned to non-existent cell labels")

    # Prepare final data structures
    cells = props_df.merge(cell_boundaries)
    cells.sort_values(by=['label', 'x_cell', 'y_cell'], inplace=True)
    spots = spots.merge(cells, how='left', on=['label'])

    return (
        cells.drop(columns=['coords']).rename(columns={
            'x_cell': 'x0',
            'y_cell': 'y0'
        }),
        cells[['label', 'coords']].rename(columns={
            'label': 'cell_id'
        }),
        spots[['x', 'y', 'label', 'Gene', 'x_cell', 'y_cell']].rename(columns={
            'Gene': 'target',
            'x': 'x_global',
            'y': 'y_global'
        })
    )
