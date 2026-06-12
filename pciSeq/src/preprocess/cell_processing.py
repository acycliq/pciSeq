""" Functions to extract the cell boundaries and calculate cell properties """
from typing import List
import numpy as np
import pandas as pd
from numba import njit
from multiprocessing import Pool, cpu_count
from multiprocessing.dummy import Pool as ThreadPool

# All that below is to avoid diplib to show a welcome msg on import. If hasattr(sys,'ps1') return True, then
# import diplib will print out the diplib name, version, description which I find annoying. I am deleting ps1
# and then reinstate it after the import.
try:
    import sys
    ps1 = sys.__dict__['ps1']
    del sys.ps1
    import diplib as dip
    sys.__dict__['ps1'] = ps1
except KeyError:
    import diplib as dip


def extract_borders_dip(label_image, offset_x=0, offset_y=0, exclude_labels=(0,)):
    """
    Extracts the cell boundaries from the label image array. The background is
    assumed to have label=0 and it will be ignored by default.
    Parameters
    ----------
    label_image:    The label image array, typically obtained from some image segmentation
                    application and maps every pixel on the image to a cell label.
    offset_x:       Amount to shift the boundaries along the x-axis
    offset_y:       Amount to shift the boundaries along the y-axis
    exclude_labels: Array-like, contains the labels to be ignored.

    Returns
    -------
    Returns a dataframe with columns ['labels', 'coords'] where column 'coords' keeps a
    list like [[x0, y0], [x1, y2],...,[x0, y0]] of the (closed-loop) boundaries coordinates
    for corresponding cell  label
    """

    if exclude_labels is None:
        exclude_labels = [0]
    labels = sorted(set(label_image.flatten()) - set(exclude_labels))
    cc = dip.GetImageChainCodes(label_image)  # input must be an unsigned integer type
    d = {}
    for c in cc:
        if c.objectID in labels:
            # p = np.array(c.Polygon())
            p = c.Polygon().Simplify()
            p = p + np.array([offset_x, offset_y])
            p = np.uint64(p).tolist()
            p.append(p[0])  # append the first pair at the end to close the polygon
            d[np.uint64(c.objectID)] = p
        else:
            pass
    df = pd.DataFrame([d]).T
    df = df.reset_index()
    df.columns = ['label', 'coords']
    return df


def extract_borders(cell_labels):
    '''
    Extracts the cell boundaries from the label image array. Same as 'extract_borders_dip()' but a lot faster.
    Parameters
    ----------
    label_image:    The label image array, typically obtained from some image segmentation
                    application and maps every pixel on the image to a cell label.

    Returns
    -------
    Returns a dataframe with columns 'labels' and 'coords'
    """
    '''
    cell_boundaries = pd.DataFrame()
    borders_list = _extract_borders(cell_labels)
    d = dict(borders_list)
    cell_boundaries['label'] = d.keys()
    cell_boundaries['coords'] = d.values()
    return cell_boundaries


def _extract_borders(label_image):
    """
    Extracts the cell boundaries from the label image array.
    Returns a dict where keys are the cell label and values the corresponding cell boundaries
    """

    # labels = sorted(set(label_image.flatten()) - set(exclude_labels))
    cc = dip.GetImageChainCodes(label_image)  # input must be an unsigned integer type

    pool = ThreadPool(cpu_count())
    # it would be nice to process only those cc whose cc.objectID is in labels
    results = pool.map(parse_chaincode, cc)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    return dict(results)


def parse_chaincode(c):
    p = c.Polygon().Simplify()
    p = np.uint64(p).tolist()
    p.append(p[0])  # append the first pair at the end to close the polygon
    return np.uint64(c.objectID), p


@njit
def _accumulate_cell_props(all_lab, all_row, all_col, all_z, n):
    """Single pass: accumulate counts, coordinate sums and z extents."""
    counts = np.zeros(n, dtype=np.int64)
    sum_x = np.zeros(n, dtype=np.float64)
    sum_y = np.zeros(n, dtype=np.float64)
    sum_z = np.zeros(n, dtype=np.float64)
    z_min = np.full(n, 999999, dtype=np.int32)
    z_max = np.full(n, -1, dtype=np.int32)

    for i in range(len(all_lab)):
        lab = all_lab[i]
        counts[lab] += 1
        sum_x[lab] += all_col[i]
        sum_y[lab] += all_row[i]
        sum_z[lab] += all_z[i]
        z_val = all_z[i]
        if z_val < z_min[lab]:
            z_min[lab] = z_val
        if z_val > z_max[lab]:
            z_max[lab] = z_val

    return counts, sum_x, sum_y, sum_z, z_min, z_max


def calculate_cell_properties(coo_list: List, voxel_size: List[float]) -> pd.DataFrame:
    """
    Calculate cell properties directly from sparse matrices.
    Gives identical results to skimage.regionprops but without
    converting to dense arrays.

    To verify against regionprops::

        masks = np.stack([coo.toarray() for coo in coo_list])
        scaling = [voxel_size[2]/voxel_size[0],
                   voxel_size[1]/voxel_size[0],
                   voxel_size[0]/voxel_size[0]]
        props = skmeas.regionprops_table(
            masks, spacing=scaling,
            properties=['label', 'area', 'centroid', 'bbox'])
        df = pd.DataFrame(props)
        # mean area per plane: divide volume by the number of planes the cell spans
        df['area'] = df['area'] / (df['bbox-3'] - df['bbox-0'])

    Parameters
    ----------
    coo_list : List[coo_matrix]
        List of sparse matrices containing cell labels, one per z-plane
    voxel_size : List[float]
        Physical size of voxels [x, y, z]

    Returns
    -------
    pd.DataFrame
        Cell properties: label, area (mean per slice), z_cell, y_cell, x_cell
    """
    scaling = [voxel_size[0] / voxel_size[0], voxel_size[1] / voxel_size[0], voxel_size[2] / voxel_size[0]]
    scaling = scaling[::-1]  # zyx order
    sz, sy, sx = scaling

    max_label = max(coo.data.max() for coo in coo_list if coo.nnz > 0)
    n = max_label + 1
    total_nnz = sum(coo.nnz for coo in coo_list)

    # Concat all sparse data into flat arrays
    all_lab = np.empty(total_nnz, dtype=np.int32)
    all_col = np.empty(total_nnz, dtype=np.float32)
    all_row = np.empty(total_nnz, dtype=np.float32)
    all_z = np.empty(total_nnz, dtype=np.int32)

    offset = 0
    for plane_idx, coo in enumerate(coo_list):
        k = coo.nnz
        if k == 0:
            continue
        s = slice(offset, offset + k)
        all_lab[s] = coo.data
        all_col[s] = coo.col
        all_row[s] = coo.row
        all_z[s] = plane_idx
        offset += k

    # Single compiled pass over all sparse data
    counts, sum_x, sum_y, sum_z, z_min, z_max = _accumulate_cell_props(
        all_lab[:offset], all_row[:offset], all_col[:offset], all_z[:offset], n
    )

    # Extract valid labels (skip background)
    valid = counts > 0
    valid[0] = False
    lab_ids = np.where(valid)[0]

    c = counts[lab_ids].astype(np.float64)
    z_extent = (z_max[lab_ids] - z_min[lab_ids] + 1).astype(np.float64)
    mean_area = (c * (sz * sy * sx)) / z_extent

    return pd.DataFrame({
        'label': lab_ids.astype(np.uint32),
        'area': mean_area.astype(np.uint32),
        'z_cell': (sum_z[lab_ids] / c * sz).astype(np.float32),
        'y_cell': (sum_y[lab_ids] / c * sy).astype(np.float32),
        'x_cell': (sum_x[lab_ids] / c * sx).astype(np.float32),
    })
