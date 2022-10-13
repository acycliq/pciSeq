"""
From https://github.com/jrussell25/dask-regionprops/tree/main/dask_regionprops
Edited to work with the label image as a list of coo matrices
"""

import re
from itertools import product

import dask.dataframe as dd
import os
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
import xarray as xr
from dask import delayed
from skimage.measure import regionprops_table
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger(__name__)

DEFAULT_PROPERTIES = (
    "label",
    "bbox",
    "centroid",
    "area",
    "convex_area",
    "eccentricity",
    "equivalent_diameter",
    "euler_number",
    "extent",
    "feret_diameter_max",
    "perimeter_crofton",
    "solidity",
    "moments_hu",
)

DEFAULT_WEIGHTED_PROPERTIES = (
    *DEFAULT_PROPERTIES,
    "centroid_weighted",
    "intensity_max",
    "intensity_mean",
    "intensity_min",
    "moments_weighted_hu",
)

DEFAULT_META = pd.DataFrame(
    regionprops_table(
        label_image=np.ones((1, 1), dtype="uint8"),
        intensity_image=np.ones((1, 1)),
        properties=DEFAULT_WEIGHTED_PROPERTIES,
    )
)

DEFAULT_PROPS_TO_COLS = {}
for prop in DEFAULT_WEIGHTED_PROPERTIES:
    col_list = []
    for c in DEFAULT_META.columns:
        stripped = re.sub("[-0-9]", "", c)
        if stripped == prop:
            col_list.append(c)
    DEFAULT_PROPS_TO_COLS[prop] = col_list


def regionprops_df(
    labels_im, intensity_im=None, props=DEFAULT_PROPERTIES, other_cols={}
):
    df = pd.DataFrame(regionprops_table(labels_im.astype(np.uint32), intensity_im, properties=props))
    for k, v in other_cols.items():
        df[k] = v
    return df


def regionprops(labels, intensity=None, properties=DEFAULT_PROPERTIES, core_dims=None):
    """
    Loop over the frames of ds and compute the regionprops for
    each labelled image in each frame.

    Parameters
    ----------
    labels : array-like of int
        Array containing labelled regions. Background is assumed to have
        value 0 and will be ignored.
    intensity : array-like or None, Default None
        Optional intensity field to compute weighted region properties from.
    properties : str, tuple[str] default "non-image"
        Properties to compute for each region. Can pass an explicit tuple
        directly to regionprops or use one of the followings shortcuts:
        "minimal", "non-image", "all". If provided an intensity image, basic
        weighted properties will also be computed by defualt.
    core_dims : tuple[int] or tuple[str] default None
        Dimensions of input arrays that correspond to spatial (xy) dimensions of each
        image. If None, it is assumed that the final two dimensions are the
        spatial dimensions.

    Returns
    -------
    regionprops_df : dask.DataFrame
        Lazily constructed dataframe containing columns for each specified
        property.
    """

    d_regionprops = delayed(regionprops_df)

    loop_sizes = _get_loop_sizes(labels)

    if intensity is not None:
        properties = DEFAULT_WEIGHTED_PROPERTIES

    meta = _get_meta(loop_sizes, properties)

    labels_arr, intensity_arr = _get_arrays(labels, intensity)

    all_props = []

    for dims in product(*(range(v) for v in loop_sizes.values())):
        other_cols = dict(zip(loop_sizes.keys(), dims))

        if intensity_arr is not None:
            frame_props = d_regionprops(
                labels_arr[dims].toarray(), intensity_arr[dims], properties, other_cols
            )
        else:
            d = dims[0]
            if np.any(labels[0].data > 0):
                frame_props = d_regionprops(
                    labels_arr[d].toarray(), None, properties, other_cols
                )
                all_props.append(frame_props)
            else:
                logger.info('Page %d of the stack is empty' % d)

    cell_props = dd.from_delayed(all_props, meta=meta)

    return cell_props


def _get_meta(loop_sizes, properties):

    meta = pd.DataFrame()
    for prop in properties:
        meta = meta.join(DEFAULT_META[DEFAULT_PROPS_TO_COLS[prop]])

    other_cols = pd.DataFrame(columns=list(loop_sizes.keys()), dtype=int)

    return meta.join(other_cols)


def _get_loop_sizes(labels):
    if np.all([isinstance(d, coo_matrix) for d in labels]) and isinstance(labels, list):
        loop_sizes = {f"dim-{0}": len(labels)}
    else:
        raise TypeError("labels should a list of coo matrices")
    return loop_sizes


def _get_pos_core_dims(core_dims, ndim):
    pos_core_dims = []
    for d in core_dims:
        if d < 0:
            pos = ndim + d
            pos_core_dims.append(pos)
        else:
            pos_core_dims.append(d)
    return tuple(pos_core_dims)


def _get_arrays(labels, intensity):

    if intensity is None:
        intensity_arr = None
    else:
        if isinstance(intensity, xr.DataArray):
            intensity_arr = intensity.data
        else:
            intensity_arr = intensity

    if isinstance(labels, xr.DataArray):
        labels_arr = labels.data
    else:
        labels_arr = labels

    return labels_arr, intensity_arr