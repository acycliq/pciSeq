import config
import h5py
import numpy as np
from scipy.sparse import load_npz
import os
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from collections import defaultdict
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()
# logger.disabled = True



def _to_csr_matrix(i, j, n):
    """Using i and j as coo-format coordinates, return csr matrix."""
    n = int(n)
    v = np.ones_like(i)
    mat = coo_matrix((v, (i, j)), shape=(n, n))
    return mat.tocsr()


def get_dir(my_config, tile_id):
    root_dir = my_config['FOV_ROOT']
    return os.path.join(root_dir, 'tile_' + str(tile_id))


def _get_connected_labels(lst):
    '''
    find which positions in the input list have repeated values
    Example:
     If lst = [0, 4, 4] then it returns [1, 2] because 4 appears twice, at positions 1 and 2 of the input list
     if lst = [0,1,2,1,3,4,2,2] then it returns [[1, 3], [2, 6, 7]]
    :param lst:
    :return:
    '''
    output = defaultdict(list)
    # Loop once over lst, store the indices of all unique elements
    for i, el in enumerate(lst):
        output[el].append(i)
    return np.array([np.array(d) for d in output.values() if len(d) > 1])


def load_mat(filepath):
    '''
    reads a Matlab mat file and returns the CellMap
    :param filepath:
    :return:
    '''
    arrays = {}
    f = h5py.File(filepath, 'r')
    for k, v in f.items():
        arrays[k] = np.array(v)

    # CAUTION: TAKE THE TRANSPOSE. MORE DETAILS ON
    # https://www.mathworks.com/matlabcentral/answers/308303-why-does-matlab-transpose-hdf5-data
    logger.info('***** Returning the transpose of the input Matlab array *******')
    return arrays['CellMap'].T

def load_mat_2(filepath):
    x = loadmat(filepath)
    return x['CellMap']

def blockfy(a, p, q):
    '''
    Divides array a into subarrays of size p-by-q
    p: block row size
    q: block column size
    '''
    # p = my_config['p']
    # q = my_config['q']
    m = a.shape[0]  # image row size (ie, dimension size left to right)
    n = a.shape[1]  # image column size (ie, dimension size top to bottom)

    # pad array with NaNs so it can be divided by p row-wise and by q column-wise
    bpr = ((n - 1) // p + 1)  # blocks per row, how many blocks left to right
    bpc = ((m - 1) // q + 1)  # blocks per column, how many blocks top to bottom
    M = p * bpr
    N = q * bpc

    A = np.nan * np.ones([N, M])
    A[:a.shape[0], :a.shape[1]] = a

    block_list = []
    previous_row = 0
    for row_block in range(bpc):
        previous_row = row_block * p
        # previous_column = 0
        for column_block in range(bpr):
            previous_column = column_block * q
            block = A[previous_row:previous_row + p, previous_column:previous_column + q]

            # remove nan columns and nan rows
            nan_cols = np.all(np.isnan(block), axis=0)
            block = block[:, ~nan_cols]
            nan_rows = np.all(np.isnan(block), axis=1)
            block = block[~nan_rows, :]

            ## append
            if block.size:
                block_list.append(block.astype(int))

    return block_list


def tilefy(a, p, q):
    """
    Splits array '''a''' into smaller subarrays of size p-by-q
    Parameters
    ----------
    a: Numpy array of size m-by-n
    p: int
    q: int

    Returns
    -------
    A list of numpy arrays

    Example:
    a = np.arange(21).reshape(7,3)
    out = tilefy(a, 4, 2)

    will return
     [array([[ 0,  1],
        [ 3,  4],
        [ 6,  7],
        [ 9, 10]]),
     array([[ 2],
            [ 5],
            [ 8],
            [11]]),
     array([[12, 13],
            [15, 16],
            [18, 19]]),
     array([[14],
            [17],
            [20]])]

    """
    m, n = a.shape
    ltr = -(-n//q)  # left to right
    ttb = -(-m//p)  # top to bottom
    out = []
    for j in range(ttb):
        rows = np.arange(p) + j*p
        for i in range(ltr):
            cols = np.arange(q) + i*q
            _slice = a[rows[0]:rows[-1]+1, cols[0]:cols[-1]+1]
            out.append(_slice.astype(np.int32))
    return out


def split_CellMap(filepath, p, q=None):
    if q is None:
        q = p
    cellmap = load_mat(filepath)
    # p = q = 2000
    out = blockfy(cellmap, p, q)
    return out


def split_label_img(filepath, p, q=None):
    """
    Splits the label_image into smaller chunks(tiles) of size p-by-q
    Parameters
    ----------
    filepath: Path to the npy file (the output of the cell segmentation)
    p: width in pixels of the tile
    q: height in pixels of the tile

    Note: keep p = q. Code has not been tested for p != q

    Returns
    -------
    """
    label_img = load_npz(filepath).toarray()
    out = tilefy(label_img, p, q)
    return out

if __name__ == "__main__":
    p = q = 2000
    filepath = os.path.join(config.ROOT_DIR, 'CellMap_left.mat')
    fov_list = split_CellMap(filepath, p, q)
    print('Done!')

