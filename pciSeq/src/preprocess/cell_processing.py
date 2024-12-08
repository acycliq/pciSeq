""" Functions to extract the cell boundaries """
from typing import List
import numpy as np
import pandas as pd
import skimage.measure as skmeas
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
    p = p + np.array([0, 0])
    p = np.uint64(p).tolist()
    p.append(p[0])  # append the first pair at the end to close the polygon
    return np.uint64(c.objectID), p


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
