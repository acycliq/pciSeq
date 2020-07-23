import config
import h5py
import numpy as np
import os
from scipy.io import loadmat
import logging
import time

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()
# logger.disabled = True


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


def split_CellMap(filepath, p, q=None):
    if q is None:
        q = p
    cellmap = load_mat(filepath)
    # p = q = 2000
    out = blockfy(cellmap, p, q)
    return out


if __name__ == "__main__":
    p = q = 2000
    filepath = os.path.join(config.ROOT_DIR, 'CellMap_left.mat')
    fov_list = split_CellMap(filepath, p, q)
    print('Done!')

