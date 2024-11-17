"""
Cell label processing functionality for pciSeq.
Handles label validation, reordering, and cell assignment.
"""

from typing import List, Union, Optional, Tuple
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
from multiprocessing.dummy import Pool as ThreadPool
import logging

label_logger = logging.getLogger(__name__)


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
        label_logger.info('Relabelling segmentation array to sequential integers...')
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
