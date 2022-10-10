import sys
import tempfile, shutil
from tqdm import tqdm
from urllib.parse import urlparse
from urllib.request import urlopen
import numpy as np
import pandas as pd
import numexpr as ne
from pciSeq.src.cell_call.log_config import logger
# import numba as nb
import os
import glob
import sqlite3
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
# logger = logging.getLogger(__name__)


# def read_image_objects(img_obj, cfg):
#     meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
#     relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius
#
#     # append 1 for the misreads
#     relCellRadius = np.append(1, relCellRadius)
#
#     nom = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(cfg['InsideCellBonus'])) + np.exp(cfg['InsideCellBonus'])
#     denom = np.exp(-0.5) * (1 - np.exp(cfg['InsideCellBonus'])) + np.exp(cfg['InsideCellBonus'])
#     CellAreaFactor = nom / denom
#
#     out = {}
#     out['area_factor'] = CellAreaFactor
#     # out['area_factor'] = np.ones(CellAreaFactor.shape)
#     # logger.info('Overriden CellAreaFactor = 1')
#     out['rel_radius'] = relCellRadius
#     out['area'] = np.append(np.nan, img_obj.area)
#     out['x'] = np.append(-sys.maxsize, img_obj.x.values)
#     out['y'] = np.append(-sys.maxsize, img_obj.y.values)
#     out['cell_label'] = np.append(0, img_obj.label.values)
#     # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
#     # and will be used to get all the misreads. It was given the label=0 and some very small
#     # negative coords
#
#     return out


def negBinLoglik(x, r, p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''

    # sanity check
    # assert (np.all(da_x.coords['cell_id'].data == da_p.coords['cell_id'])), 'gene counts and beta probabilities are not aligned'
    # assert (np.all(da_x.coords['gene_name'].data == da_p.coords['gene_name'])), 'gene counts and beta probabilities are not aligned'

    contr = np.zeros(p.shape)
    x = x[:, :, None]
    ne.evaluate("x * log(p) + r * log(1 - p)", out=contr)
    return contr


# @nb.njit(parallel=True, fastmath=True)
# def nb_negBinLoglik(x, r, p):
#     '''
#     Negative Binomial loglikehood
#     :param x:
#     :param r:
#     :param p:
#     :return:
#     '''
#     out = np.empty(p.shape,p.dtype)
#
#     for i in nb.prange(p.shape[0]):
#         for j in range(p.shape[1]):
#             if x[i, j, 0] != 0.:
#                 x_ = x[i, j, 0]
#                 for k in range(p.shape[2]):
#                     out[i, j, k] = x_ * np.log(p[i, j, k]) + r * np.log(1.-p[i, j, k])
#             else:
#                 for k in range(p.shape[2]):
#                     out[i, j, k] = r * np.log(1.-p[i, j, k])
#
#     return out



def softmax(X, theta = 1.0, axis = None):
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
    y = y - np.expand_dims(np.max(y, axis = axis), axis)

    # exponentiate y
    y = np.exp(y)

    # take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis = axis), axis)

    # finally: divide elementwise
    p = y / ax_sum

    # flatten if X was 1D
    if len(X.shape) == 1: p = p.flatten()

    return p


def hasConverged(spots, p0, tol):
    p1 = spots.parent_cell_prob
    if p0 is None:
        p0 = np.zeros(p1.shape)
    delta = np.max(np.abs(p1 - p0))
    converged = (delta < tol)
    return converged, delta


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
        if size > mb_size*1024*1024:
            logger.info('saved %s with file size %4.3f MB' % (file_out, size/(1024*1024)))
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
        files = glob.glob(OUT_DIR + '/*.'+ext)
        for f in files:
            os.remove(f)

    for i, d in enumerate(df_list):
        fname = os.path.join(OUT_DIR, filename + '_%d.%s' % (i, ext))
        if ext == 'json':
            d.to_json(fname,  orient='records')
        elif ext == 'tsv':
            d.to_csv(fname, sep='\t', index=False)


def _get_file(OUT_DIR, filepath, n, header_line):
    [filename, ext] = os.path.basename(filepath).split('.')
    file = os.path.join(OUT_DIR, filename + '_%d.%s' % (n, ext))
    handle = open(file, "a")
    handle.write(header_line)
    return file, handle


def download_url_to_file(url, dst, progress=True):
    r"""Download object at the given URL to a local path.
            Thanks to torch, slightly modified
    Args:
        url (string): URL of the object to download
        dst (string): Full path where object will be saved, e.g. `/tmp/temporary_file`
        progress (bool, optional): whether or not to display a progress bar to stderr
            Default: True
    """
    file_size = None
    u = urlopen(url)
    meta = u.info()
    if hasattr(meta, 'getheaders'):
        content_length = meta.getheaders("Content-Length")
    else:
        content_length = meta.get_all("Content-Length")
    if content_length is not None and len(content_length) > 0:
        file_size = int(content_length[0])
    # We deliberately save it in a temp file and move it after
    dst = os.path.expanduser(dst)
    dst_dir = os.path.dirname(dst)
    f = tempfile.NamedTemporaryFile(delete=False, dir=dst_dir)
    try:
        with tqdm(total=file_size, disable=not progress,
                  unit='B', unit_scale=True, unit_divisor=1024) as pbar:
            while True:
                buffer = u.read(8192)
                if len(buffer) == 0:
                    break
                f.write(buffer)
                pbar.update(len(buffer))
        f.close()
        shutil.move(f.name, dst)
    finally:
        f.close()
        if os.path.exists(f.name):
            os.remove(f.name)


def load_from_url(url):
    # files = []
    parts = urlparse(url)
    filename = os.path.basename(parts.path)
    if not os.path.exists(filename):
        sys.stderr.write('Downloading: "{}" to {}\n'.format(url, filename))
        download_url_to_file(url, filename)
    return filename


def gaussian_ellipsoid(mu, cov, sdwidth=None):
    if sdwidth is None:
        sdwidth = 1

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

    return np.array(list(zip(*out)))


def gaussian_ellipsoid_props(cov, sdwidth=None):
    """
    get the scaling, rotation of the ellipsoid
    """
    tol = 1.0e-10
    cov = np.where(cov < tol, 0, cov)
    eigvals, eigvecs = np.linalg.eig(cov)
    scaling = sdwidth * np.sqrt(eigvals)
    # rotation = roll_pitch_yaw(eigvecs)
    rotation = euler_angles(eigvecs.T)
    return scaling, rotation


def euler_angles(r):
    theta_x = np.arctan2(r[2, 1], r[2, 2])
    theta_y = np.arcsin(r[2,0])
    theta_z = np.arctan2(r[1, 0], r[0, 0])

    return [theta_x, theta_y, theta_z]


def get_db_tables(con):
    str = "SELECT name FROM sqlite_schema WHERE type = 'table' ORDER BY name;"
    tables = con.execute(str).fetchall()
    if len(tables) > 0:
        return [d[0] for d in tables]
    else:
        return tables


def drop_db_table(table_name, con):
    cmd = "DROP TABLE %s" % table_name
    con.execute(cmd).fetchall()
    logger.info('Table %s dropped from the database' % table_name)


def db_connect(dbpath, remove_if_exists=True):
    if os.path.isfile(dbpath) and remove_if_exists:
        try:
            os.remove(dbpath)
        except PermissionError:
            con = sqlite3.connect(dbpath, uri = True)
            tables = get_db_tables(con)
            for t in tables:
                drop_db_table(t, con)
            return con

    return sqlite3.connect(dbpath, uri = True)


def get_out_dir(path, sub_folder=''):
    if path[0] == 'default':
        out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq', sub_folder)
    else:
        out_dir = os.path.join(path[0], sub_folder)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


