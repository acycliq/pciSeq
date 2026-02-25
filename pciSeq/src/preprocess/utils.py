"""
Logging utilities for preprocessing operations.
"""

from typing import List, Tuple
import logging
import numpy as np
from scipy.sparse import coo_matrix
import pandas as pd

logger = logging.getLogger(__name__)


def log_data_summary(spots: pd.DataFrame,
                     coo: List[coo_matrix],
                     dimensions: Tuple[int, int, int]) -> None:
    """
    Log summary of data dimensions and counts.

    Args:
        spots: DataFrame with spot data
        coo: List of sparse matrices containing cell labels
        dimensions: (n_planes, height, width) of image
    """
    n, h, w = dimensions

    logger.info(f'Number of spots passed-in: {spots.shape[0]}')
    logger.info(
        f'Number of segmented cells: {max([d.data.max() for d in coo if len(d.data) > 0])}')

    if n == 1:
        logger.info(f'Image dimensions: {w}px × {h}px')
    else:
        logger.info(f'Image dimensions: {n} planes, {w}px × {h}px')
