"""Utilities for processing and manipulating cell data."""
import sys
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any, Optional, Union
from numbers import Number
import logging

# Configure logging
cell_utils_logger = logging.getLogger(__name__)


def read_image_objects(img_obj: pd.DataFrame,
                       cfg: Dict[str, Any]) -> Tuple[Dict[str, np.ndarray], np.float32]:
    """Process image objects and calculate cell properties.

    Args:
        img_obj: DataFrame containing cell image objects
        cfg: Configuration dictionary

    Returns:
        Tuple containing:
            - Dict of cell properties (area_factor, rel_radius, area, x0, y0, z0, cell_label)
            - Mean cell radius as float32

    Notes:
        Handles special case of misreads by appending dummy cell with label=0
    """
    # Calculate mean cell radius and relative radii
    meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
    relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

    # Append 1 for the misreads
    relCellRadius = np.append(1, relCellRadius)

    # Calculate cell area factor
    InsideCellBonus = cfg['InsideCellBonus'] if cfg['InsideCellBonus'] else 0

    numer = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    denom = np.exp(-0.5) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    CellAreaFactor = numer / denom

    # Build output dictionary
    out = {
        'area_factor': CellAreaFactor.astype(np.float32),
        'rel_radius': relCellRadius.astype(np.float32),
        'area': np.append(np.nan, img_obj.area.astype(np.uint32)),
        'x0': np.append(-sys.maxsize, img_obj.x0.values).astype(np.float32),
        'y0': np.append(-sys.maxsize, img_obj.y0.values).astype(np.float32),
        'z0': np.append(-sys.maxsize, img_obj.z0.values).astype(np.float32),
        'cell_label': np.append(0, img_obj.label.values).astype(np.uint32)
    }

    # Add old labels if present
    if 'old_label' in img_obj.columns:
        out['cell_label_old'] = np.append(0, img_obj.old_label.values).astype(np.uint32)

    return out, meanCellRadius.astype(np.float32)


def recover_original_labels(cellData: pd.DataFrame,
                            geneData: pd.DataFrame,
                            label_map: Optional[Dict[int, int]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Restore original cell labels using label mapping.

    Args:
        cellData: Cell data DataFrame
        geneData: Gene data DataFrame
        label_map: Dictionary mapping new labels to original labels

    Returns:
        Tuple of (updated cellData, updated geneData)
    """
    if label_map is None:
        return cellData, geneData

    # Create reverse mapping
    reverse_map = {v: k for k, v in label_map.items()}

    # Update cell numbers
    cellData = cellData.assign(
        Cell_Num=cellData.Cell_Num.map(lambda x: reverse_map.get(x, x))
    )

    # Update gene data neighbors
    geneData = geneData.assign(
        neighbour=geneData.neighbour.map(lambda x: fetch_label(x, reverse_map)),
        neighbour_array=geneData.neighbour_array.map(lambda x: fetch_label(x, reverse_map))
    )
    cell_utils_logger.info("Restored original cell segmentation labels")
    return cellData, geneData


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

