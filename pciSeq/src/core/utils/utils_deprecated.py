from typing import Tuple, Optional, Dict, Any, Union
import numbers
import os
import glob
import pickle
import sys
import tempfile
from pathlib import Path
import shutil
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from tqdm import tqdm

from pciSeq import config
from pciSeq.src.diagnostics.utils import check_redis_server
from pciSeq.src.viewer.utils import (copy_viewer_code, make_config_js,
                                     make_classConfig_js, make_glyphConfig_js,
                                     make_classConfig_nsc_js, build_pointcloud,
                                     cellData_rgb, cell_gene_counts)
import logging


utils_logger = logging.getLogger(__name__)


# def init(opts):
#     """
#     Reads the opts dict and if not None, it will override the default parameter value by
#     the value that the dictionary key points to.
#     If opts is None, then the default values as these specified in the config.py file
#     are used without any change.
#     """
#     cfg = config.DEFAULT
#     log_file(cfg)
#     cfg['is_redis_running'] = check_redis_server()
#     if opts is not None:
#         default_items = set(cfg.keys())
#         user_items = set(opts.keys())
#         assert user_items.issubset(default_items), ('Options passed-in should be a dict with keys: %s ' % default_items)
#         for item in opts.items():
#             if isinstance(item[1], (int, float, list, str, dict)) or isinstance(item[1](1), np.floating):
#                 val = item[1]
#             # elif isinstance(item[1], list):
#             #     val = item[1]
#             else:
#                 raise TypeError("Only integers, floats and lists are allowed")
#             cfg[item[0]] = val
#             utils_logger.info('%s is set to %s' % (item[0], cfg[item[0]]))
#     return cfg
#
#
# def log_file(cfg):
#     """
#     Setup the logger file handler if it doesn't already exist.
#     """
#     root_logger = logging.getLogger()
#     if root_logger.handlers:
#         # setup a FileHandler if it has not been setup already. Maybe I should be adding a FileHandler anyway,
#         # regardless whether there is one already or not
#         if not np.any([isinstance(d, logging.FileHandler) for d in root_logger.handlers]):
#             logfile = os.path.join(get_out_dir(cfg['output_path']), 'pciSeq.log')
#             fh = logging.FileHandler(logfile, mode='w')
#             formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#             fh.setFormatter(formatter)
#
#             root_logger.addHandler(fh)
#             utils_logger.info('Writing to %s' % logfile)


def negative_binomial_loglikelihood(x: np.ndarray, r: float, p: np.ndarray) -> np.ndarray:
    """
    Calculate the Negative Binomial log-likelihood for given parameters.

    The Negative Binomial distribution models the number of successes (x) before
    r failures occur, with probability of success p. The PMF is:
    P(X=k) = C(k+r-1,k) * p^k * (1-p)^r

    Args:
        x: Array of observed counts with shape
        r: Number of failures until stopping (dispersion parameter)
        p: Probability of success, array with shape

    Returns:
        Array of log-likelihood contributions with shape
    """
    try:
        x = x[:, :, None]  # Add dimension for broadcasting

        # negative binomial log-likelihood.
        # Focusing only on the terms that involve p and r (without the binomial coefficient):
        log_likelihood = x * np.log(p) + r * np.log(1 - p)

        return log_likelihood

    except Exception as e:
        utils_logger.error(f"Error calculating negative binomial log-likelihood: {str(e)}")
        raise ValueError("Failed to compute log-likelihood. Check input dimensions and values.")


def softmax(X, theta=1.0, axis=None):
    """
    From https://nolanbconaway.github.io/blog/2017/softmax-numpy
    Compute the softmax of each element along an axis of X.

    Parameters
    ----------
    X: ND-Array. Probably should be floats.
    theta (optional): float parameter, used as a multiplier
        prior to exponentiation. Default = 1.0
    axis (optional): axis to compute values along. Default is the
        first non-singleton axis.

    Returns an array the same size as X. The result will sum to 1
    along the specified axis.
    """

    # make X at least 2d
    y = np.atleast_2d(X)

    # find axis
    if axis is None:
        axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)

    # multiply y against the theta parameter,
    y = y * float(theta)

    # subtract the max for numerical stability
    y = y - np.expand_dims(np.max(y, axis=axis), axis)

    # exponentiate y
    y = np.exp(y)

    # take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis=axis), axis)

    # finally: divide elementwise
    p = y / ax_sum

    # flatten if X was 1D
    if len(X.shape) == 1: p = p.flatten()

    return p


def has_converged(
        spots: Any,
        p0: Optional[np.ndarray],
        tol: float
) -> Tuple[bool, float]:
    """
    Check if the probability assignments have converged.

    Args:
        spots: Spot data object containing parent_cell_prob
        p0: Previous probability matrix (None for first iteration)
        tol: Convergence tolerance threshold

    Returns:
        Tuple containing:
            - bool: True if converged, False otherwise
            - float: Maximum absolute difference between iterations
    """
    p1 = spots.parent_cell_prob
    if p0 is None:
        p0 = np.zeros_like(p1)

    try:
        delta = np.max(np.abs(p1 - p0))
        converged = (delta < tol)
        return converged, delta
    except Exception as e:
        utils_logger.error(f"Convergence check failed: {str(e)}")
        raise


def splitter_mb(filepath, mb_size):
    """ Splits a text file in (almost) equally sized parts on the disk. Assumes that there is a header in the first line
    :param filepath: The path of the text file to be broken up into smaller files
    :param mb_size: size in MB of each chunk
    :return:
    """
    handle = open(filepath, 'r')
    OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.*')
        for f in files:
            os.remove(f)

    n = 0
    size = None
    header_line = next(handle)
    file_out, handle_out = _get_file(OUT_DIR, filepath, n, header_line)
    for line in handle:
        size = os.stat(file_out).st_size
        if size > mb_size * 1024 * 1024:
            utils_logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
            n += 1
            handle_out.close()
            file_out, handle_out = _get_file(OUT_DIR, filepath, n, header_line)
        handle_out.write(str(line))

    # print(str(file_out) + " file size = \t" + str(size))
    print('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
    handle_out.close()


def splitter_n(filepath, n):
    """ Splits a text file in n smaller files
    :param filepath: The path of the text file to be broken up into smaller files
    :param n: determines how many smaller files will be created
    :return:
    """
    filename_ext = os.path.basename(filepath)
    [filename, ext] = filename_ext.split('.')

    OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if ext == 'json':
        df = pd.read_json(filepath)
    elif ext == 'tsv':
        df = pd.read_csv(filepath, sep='\t')
    else:
        df = None

    df_list = np.array_split(df, n)
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.' + ext)
        for f in files:
            os.remove(f)

    for i, d in enumerate(df_list):
        fname = os.path.join(OUT_DIR, filename + '_%d.%s' % (i, ext))
        if ext == 'json':
            d.to_json(fname, orient='records')
        elif ext == 'tsv':
            d.to_csv(fname, sep='\t', index=False)


def _get_file(OUT_DIR, filepath, n, header_line):
    [filename, ext] = os.path.basename(filepath).split('.')
    file = os.path.join(OUT_DIR, filename + '_%d.%s' % (n, ext))
    handle = open(file, "a")
    handle.write(header_line)
    return file, handle


# def download_url_to_file(url, dst, progress=True):
#     r"""Download object at the given URL to a local path.
#             Thanks to torch, slightly modified
#     Args:
#         url (string): URL of the object to download
#         dst (string): Full path where object will be saved, e.g. `/tmp/temporary_file`
#         progress (bool, optional): whether or not to display a progress bar to stderr
#             Default: True
#     """
#     file_size = None
#     u = urlopen(url)
#     meta = u.info()
#     if hasattr(meta, 'getheaders'):
#         content_length = meta.getheaders("Content-Length")
#     else:
#         content_length = meta.get_all("Content-Length")
#     if content_length is not None and len(content_length) > 0:
#         file_size = int(content_length[0])
#     # We deliberately save it in a temp file and move it after
#     dst = os.path.expanduser(dst)
#     dst_dir = os.path.dirname(dst)
#     f = tempfile.NamedTemporaryFile(delete=False, dir=dst_dir)
#     try:
#         with tqdm(total=file_size, disable=not progress,
#                   unit='B', unit_scale=True, unit_divisor=1024) as pbar:
#             while True:
#                 buffer = u.read(8192)
#                 if len(buffer) == 0:
#                     break
#                 f.write(buffer)
#                 pbar.update(len(buffer))
#         f.close()
#         shutil.move(f.name, dst)
#     finally:
#         f.close()
#         if os.path.exists(f.name):
#             os.remove(f.name)
#
#
# def load_from_url(url):
#     # files = []
#     parts = urlparse(url)
#     filename = os.path.basename(parts.path)
#     if not os.path.exists(filename):
#         sys.stderr.write('Downloading: "{}" to {}\n'.format(url, filename))
#         download_url_to_file(url, filename)
#     return filename


# def get_out_dir(path: Optional[str] = None, sub_folder: str = '') -> str:
#     """
#     Get or create output directory path.
#
#     Args:
#         path: List containing base path, or None for default temp directory
#         sub_folder: Optional subdirectory name
#
#     Returns:
#         str: Path to output directory
#
#     Notes:
#         - Uses system temp directory if path is None or ['default']
#         - Creates directories if they don't exist
#     """
#     if path is None or path == 'default':
#         out_dir = Path(tempfile.gettempdir()) / 'pciSeq'
#     else:
#         out_dir = Path(path) / sub_folder / 'pciSeq'
#
#     out_dir.mkdir(parents=True, exist_ok=True)
#     return str(out_dir)


def gaussian_contour(mu, cov, sdwidth=3):
    """
    Draws an ellipsoid for a given covariance matrix cov
    and mean vector mu

    Example
    cov_1 = [[1, 0.5], [0.5, 1]]
    means_1 = [1, 1]

    cov_2 = [[1, -0.7], [-0.7, 1]]
    means_2 = [2, 1.5]

    cov_3 = [[1, 0], [0, 1]]
    means_3 = [0, 0]

    ellipsis_1 = gaussian_contour(means_1, cov_1)
    ellipsis_2 = gaussian_contour(means_2, cov_2)
    ellipsis_3 = gaussian_contour(means_3, cov_3)
    ellipsis_3b = gaussian_contour(means_3, cov_3, sdwidth=2)
    ellipsis_3c = gaussian_contour(means_3, cov_3, sdwidth=3)

    plt.plot(ellipsis_1[0], ellipsis_1[1], c='b')
    plt.plot(ellipsis_2[0], ellipsis_2[1], c='r')
    plt.plot(ellipsis_3[0], ellipsis_3[1], c='g')
    plt.plot(ellipsis_3b[0], ellipsis_3b[1], c='g')
    plt.plot(ellipsis_3c[0], ellipsis_3c[1], c='g')
    """

    # cov_00 = sigma_x * sigma_x
    # cov_10 = rho * sigma_x * sigma_y
    # cov_11 = sigma_y * sigma_y
    # cov = np.array([[cov_00, cov_10], [cov_10, cov_11]])
    mu = np.array(mu)

    npts = 40
    tt = np.linspace(0, 2 * np.pi, npts)
    ap = np.zeros((2, npts))
    x = np.cos(tt)
    y = np.sin(tt)
    ap[0, :] = x
    ap[1, :] = y

    eigvals, eigvecs = np.linalg.eig(cov)
    eigvals = sdwidth * np.sqrt(eigvals)
    eigvals = eigvals * np.eye(2)

    vd = eigvecs.dot(eigvals)
    out = vd.dot(ap) + mu.reshape(2, -1)

    return np.array(list(zip(*out)), dtype=np.float32)


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

    out = {}
    out['area_factor'] = CellAreaFactor.astype(np.float32)
    out['rel_radius'] = relCellRadius.astype(np.float32)
    out['area'] = np.append(np.nan, img_obj.area.astype(np.uint32))
    out['x0'] = np.append(-sys.maxsize, img_obj.x0.values).astype(np.float32)
    out['y0'] = np.append(-sys.maxsize, img_obj.y0.values).astype(np.float32)
    out['z0'] = np.append(-sys.maxsize, img_obj.z0.values).astype(np.float32)
    out['cell_label'] = np.append(0, img_obj.label.values).astype(np.uint32)
    if 'old_label' in img_obj.columns:
        out['cell_label_old'] = np.append(0, img_obj.old_label.values).astype(np.uint32)
    # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
    # and will be used to get all the misreads. It was given the label=0 and some very small
    # negative coords

    return out, meanCellRadius.astype(np.float32)


def keep_labels_unique(scdata):
    """
    In the single cell data you might find cases where two or more rows have the same gene label
    In these cases keep the row with the highest total gene count
    """

    # 1. get the row total and assign it to a new column
    scdata = scdata.assign(total=scdata.sum(axis=1))

    # 2. rank by gene label and total gene count and keep the one with the highest total
    scdata = scdata.sort_values(['gene_name', 'total'], ascending=[True, False]).groupby('gene_name').head(1)

    # 3. Drop the total column and return
    return scdata.drop(['total'], axis=1)


# @dask.delayed
def scaled_exp(cell_area_factor, sc_mean_expressions, inefficiency):
    if np.all(cell_area_factor == 1):
        subscripts = 'gk, g -> gk'
        operands = [sc_mean_expressions, inefficiency]
    else:
        subscripts = 'c, gk, g -> cgk'
        operands = [cell_area_factor, sc_mean_expressions, inefficiency]
    out = np.einsum(subscripts, *operands)
    return out


def adjust_for_anisotropy(spots, voxel_size):
    gene_col = spots.gene_name.values[:, None]
    z_plane = np.float32(spots.z_plane.values[:, None])

    data = spots[['x', 'y', 'z_plane']]
    spots_adj = anisotropy_calc(data, voxel_size)

    spots_adj = np.hstack([gene_col, spots_adj, z_plane])
    spots_adj = pd.DataFrame(spots_adj, columns=['gene_name', 'x', 'y', 'z', 'z_plane'])
    spots_adj = spots_adj.astype({'x': 'float32',
                                  'y': 'float32',
                                  'z': 'float32',
                                  'z_plane': 'float32',
                                  'gene_name': str})
    return spots_adj  # Nspots x 3


def anisotropy_calc(data, voxel_size):
    x, y, z = voxel_size
    Sx = x / x
    Sy = y / x
    Sz = z / x
    scaling_matrix = np.array([
        [Sx, 0, 0],
        [0, Sy, 0],
        [0, 0, Sz]
    ], dtype=np.float32)
    out = scaling_matrix.dot(data.T)  # 3 x Nspots
    return out.T  # Nspots x 3


def truncate_data(spots, label_image, z_min, z_max):
    if z_min is not None and z_max is not None:
        utils_logger.info(' Truncating masks and spots. Keeping those between frame %d and %d only' % (z_min, z_max))
        spots = truncate_spots(spots, z_min, z_max)
    return spots, label_image


def truncate_spots(spots, zmin, zmax):
    spots_min = spots[(spots.z_stack <= zmax) & (spots.z_stack >= zmin)]
    spots_min = spots_min.assign(z_stack=spots_min.z_stack - zmin)
    # out = spots_min.z - zmin
    return spots_min.reset_index(drop=True)


def get_img_shape(coo):
    n = len(coo)
    img_shape = set([d.shape for d in coo])
    assert len(img_shape) == 1, 'pages do not have the same shape'
    img_shape = img_shape.pop()
    w = img_shape[1]
    h = img_shape[0]
    return [n, h, w]


def gaussian_ellipsoid_props(cov, sdwidth=3):
    """
    get the scaling, rotation of the ellipsoid
    """
    tol = 1.0e-10
    cov = np.where(cov < tol, 0, cov)
    eigvals, eigvecs = np.linalg.eig(cov)
    scaling = sdwidth * np.sqrt(eigvals)
    # rotation = roll_pitch_yaw(eigvecs)
    rotation = euler_angles(eigvecs.T)
    scaling = scaling.tolist()
    return scaling, rotation


def euler_angles(r):
    theta_x = np.arctan2(r[2, 1], r[2, 2])
    theta_y = np.arcsin(r[2, 0])
    theta_z = np.arctan2(r[1, 0], r[0, 0])

    return [theta_x, theta_y, theta_z]


def gaussian_contour(mu, cov, sdwidth=3):
    """
    Draws an ellipsoid for a given covariance matrix cov
    and mean vector mu

    Example
    cov_1 = [[1, 0.5], [0.5, 1]]
    means_1 = [1, 1]

    cov_2 = [[1, -0.7], [-0.7, 1]]
    means_2 = [2, 1.5]

    cov_3 = [[1, 0], [0, 1]]
    means_3 = [0, 0]

    ellipsis_1 = gaussian_contour(means_1, cov_1)
    ellipsis_2 = gaussian_contour(means_2, cov_2)
    ellipsis_3 = gaussian_contour(means_3, cov_3)
    ellipsis_3b = gaussian_contour(means_3, cov_3, sdwidth=2)
    ellipsis_3c = gaussian_contour(means_3, cov_3, sdwidth=3)

    plt.plot(ellipsis_1[0], ellipsis_1[1], c='b')
    plt.plot(ellipsis_2[0], ellipsis_2[1], c='r')
    plt.plot(ellipsis_3[0], ellipsis_3[1], c='g')
    plt.plot(ellipsis_3b[0], ellipsis_3b[1], c='g')
    plt.plot(ellipsis_3c[0], ellipsis_3c[1], c='g')
    """

    # cov_00 = sigma_x * sigma_x
    # cov_10 = rho * sigma_x * sigma_y
    # cov_11 = sigma_y * sigma_y
    # cov = np.array([[cov_00, cov_10], [cov_10, cov_11]])
    mu = np.array(mu)

    npts = 40
    tt = np.linspace(0, 2 * np.pi, npts)
    ap = np.zeros((2, npts))
    x = np.cos(tt)
    y = np.sin(tt)
    ap[0, :] = x
    ap[1, :] = y

    eigvals, eigvecs = np.linalg.eig(cov)
    eigvals = sdwidth * np.sqrt(eigvals)
    eigvals = eigvals * np.eye(2)

    vd = eigvecs.dot(eigvals)
    out = vd.dot(ap) + mu.reshape(2, -1)

    return np.array(list(zip(*out)), dtype=np.float32)


def recover_original_labels(cellData, geneData, label_map):
    """Restore original cell labels using CellLabelManager mapping"""
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

    return cellData, geneData


def fetch_label(x, d):
    x = [x] if isinstance(x, numbers.Number) else x
    out = [d[v] for v in x]
    return out[0] if len(out) == 1 else out


# def serialise(varBayes, debug_dir):
#     if not os.path.exists(debug_dir):
#         os.makedirs(debug_dir)
#     pickle_dst = os.path.join(debug_dir, 'pciSeq.pickle')
#     with open(pickle_dst, 'wb') as outf:
#         pickle.dump(varBayes, outf)
#         utils_logger.info('Saved at %s' % pickle_dst)


def purge_spots(spots, sc):
    drop_spots = list(set(spots.gene_name) - set(sc.index))
    utils_logger.warning('Found %d genes that are not included in the single cell data' % len(drop_spots))
    idx = ~ np.in1d(spots.gene_name, drop_spots)
    spots = spots.iloc[idx]
    utils_logger.warning('Removed from spots: %s' % drop_spots)
    return spots


def export_db_tables(out_dir, con):
    tables = con.get_db_tables()
    for table in tables:
        export_db_table(table, out_dir, con)


def export_db_table(table_name, out_dir, con):
    df = con.from_redis(table_name)
    fname = os.path.join(out_dir, table_name + '.csv')
    df.to_csv(fname, index=False)
    utils_logger.info('Saved at %s' % fname)


# def write_data(cellData, geneData, cellBoundaries, varBayes, cfg):
#     dst = get_out_dir(cfg['output_path'])
#     out_dir = os.path.join(dst, 'data')
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
#
#     cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
#     utils_logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))
#
#     geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
#     utils_logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))
#
#     # Do not change the if-then branching flow. InsideCellBonus can take the value of zero
#     # and the logic below will ensure that in that case, the segmentation borders will be
#     # drawn instead of the gaussian contours.
#     if cfg['InsideCellBonus'] is False:
#         ellipsoidBoundaries = cellData[['Cell_Num', 'gaussian_contour']]
#         ellipsoidBoundaries = ellipsoidBoundaries.rename(columns={"Cell_Num": "cell_id", "gaussian_contour": "coords"})
#         ellipsoidBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
#         utils_logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellBoundaries.tsv')))
#     else:
#         cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
#         utils_logger.info('Saved at %s' % (os.path.join(out_dir, 'cellBoundaries.tsv')))
#
#     serialise(varBayes, os.path.join(out_dir, 'debug'))


def validate(spots, coo, scData, cfg):

    # check the label image
    coo = _validate_coo(coo)

    # Check now the config dict. It needs to be done after we have checked the coo
    # because the shape of the label image determines whether the data is 3D or 2D.
    cfg = _validate_cfg(cfg, coo)

    # check the spots and see if all genes are in your single cell library
    spots = _validate_spots(spots, coo, scData, cfg)

    # check single cell data
    if scData is not None:
        assert isinstance(scData,
                          pd.DataFrame), "Single cell data should be passed-in to the fit() method as a dataframe"

    # remove genes that cannot been found in the single cell data
    spots = clean_spots(spots, scData)

    # do some datatype casting
    spots = spots.astype({
        'gene_name': str,
        'x': np.float32,
        'y': np.float32,
        'z_plane': np.float32})

    if isinstance(cfg['MisreadDensity'], numbers.Number):
        val = cfg['MisreadDensity']
        cfg['MisreadDensity'] = {'default': val}
    elif isinstance(cfg['MisreadDensity'], dict):
        if 'default' not in cfg['MisreadDensity']:
            raise ValueError("When MisreadDensity is a dictionary, it must contain a 'default' key")
    else:
        raise ValueError("MisreadDensity must be either a number or a dictionary with a 'default' key")

    return spots, coo, cfg


def _validate_coo(coo):
    if isinstance(coo, coo_matrix):
        coo = [coo]
    if isinstance(coo, list):
        assert np.all([isinstance(d, coo_matrix) for d in coo]), (
            'The segmentation masks should be passed-in as a coo_matrix')

        # check all planes have the same shape
        shapes = [d.shape for d in coo]
        assert shapes.count(shapes[0]) == len(shapes), 'Not all elements in coo have the same shape'
    else:
        raise ValueError('The segmentation masks should be passed-in as a coo_matrix')
    return coo


def _validate_spots(spots, coo, sc, cfg):
    assert (isinstance(spots, pd.DataFrame)), 'Spots should be provided as a dataframe'

    # check the spots
    if not cfg['is3D']:
        assert {'gene_name', 'x', 'y'}.issubset(spots.columns), ("Spots should be passed-in to the fit() "
                                                            "method as a dataframe with columns ['gene_name', 'x', 'y']")
    else:
        assert set(spots.columns) == {'gene_name', 'x', 'y', 'z_plane'}, ("Spots should be passed-in to the fit() method "
                                                                     "as a dataframe with columns 'gene_name', 'x',  'y', "
                                                                     "'z_plane'")

    if len(coo) > 1 and set(spots.columns) == {'gene_name', 'x', 'y'}:
        raise ValueError("Label image is 3D, the spots are 2D")

    # if you have feed a 2d spot, append a flat z dimension so that everything in the code from now on will
    # be treated as 3d
    if isinstance(spots, pd.DataFrame) and set(spots.columns) == {'gene_name', 'x', 'y'}:
        spots = spots.assign(z_plane=0)

    return spots


def _validate_cfg(cfg, coo):
    # The label image will determine whether dataset is 3d or 2d
    cfg['is3D'] = False if len(coo) == 1 else True

    # override the InsideCellBonus setting and set it to false if the dataset is 3D
    # maybe drop a logger warning msg here too!
    if cfg['is3D']:
        cfg['InsideCellBonus'] = False

    # if InsideCellBonus is True, set it to a numerical value
    if cfg['InsideCellBonus'] is True:
        """
        This is not good enough! The default value for InsideCellBonus is now kept in two places, config.py and 
        here. What happens if I change the config.py and I set InsideCellBonus = 3 for example? 
        The line below will stll set it 2 which is not the default anymore! 
        """
        d = 2
        cfg['InsideCellBonus'] = d
        utils_logger.warning('InsideCellBonus was passed-in as True. Overriding with the default value of %d' % d)

    if cfg['cell_type_prior'].lower() not in ['uniform'.lower(), 'weighted'.lower()]:
        raise ValueError("'cel_type_prior' should be either 'uniform' or 'weighted' ")

    # make sure the string is lowercase from now on
    cfg['cell_type_prior'] = cfg['cell_type_prior'].lower()

    if not isinstance(cfg['launch_viewer'], bool):
        assert cfg['launch_viewer'].lower() == '2d', "'launch_viewer' should be True, False or '2d' "
        cfg['launch_viewer'] = cfg['launch_viewer'].lower()

    return cfg


def clean_spots(spots, sc):
    # remove genes that cannot been found in the single cell data
    if (sc is not None) and (not set(spots.gene_name).issubset(sc.index)):
        spots = purge_spots(spots, sc)
    return spots


def pre_launch(cellData, geneData, coo, scRNAseq, cfg):
    '''
    Does some housekeeping, pre-flight control checking before the
    viewer is triggered

    Returns the destination folder that keeps the viewer code and
    will be served to launch the website
    '''
    cfg['is3D'] = False if cfg['launch_viewer'] == '2d' else cfg['launch_viewer']
    [n, h, w] = get_img_shape(coo)
    dst = get_out_dir(cfg['output_path'])
    pciSeq_dir = copy_viewer_code(cfg, dst)
    make_config(dst, pciSeq_dir, (cellData, scRNAseq, h, w))
    if cfg['is3D']:
        # ***************************************
        # OK FOR NOW BUT THIS IS A TEMP FIX ONLY
        # THESE LINES MUST BE REMOVED
        cellData.X = cellData.X - w / 2
        cellData.Y = cellData.Y - h / 2
        cellData.Z = cellData.Z - n / 2

        geneData.x = geneData.x - w / 2
        geneData.y = geneData.y - h / 2
        geneData.z = geneData.z - n / 2
        # ***************************************

        build_pointcloud(geneData, pciSeq_dir, dst)
        cellData_rgb(cellData, pciSeq_dir, dst)
        cell_gene_counts(geneData, dst)
    return dst


def make_config(target_dir, source_dir, data):
    '''
    makes the three javascript configuration files that drive the viewer.
    They are getting saved into the dst folder. By default this is the tmp directory
    '''

    cellData = data[0]
    scRNAseq = data[1]
    img_shape = data[2:]
    # 1. make the main configuration file
    make_config_js(target_dir, img_shape)

    # 2. make the file for the class colors
    if scRNAseq is None:
        label_list = [d[:] for d in cellData.ClassName.values]
        labels = [item for sublist in label_list for item in sublist]
        labels = sorted(set(labels))

        make_classConfig_nsc_js(labels, target_dir)
    else:
        make_classConfig_js(source_dir, target_dir)

    # 3. make the file for the glyph colors
    gene_panel = scRNAseq.index.values
    make_glyphConfig_js(gene_panel, source_dir, target_dir)


