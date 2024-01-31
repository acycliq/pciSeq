from pciSeq.src.core.logger import logger_setup
from pciSeq import fit
from scipy.sparse import load_npz, coo_matrix
import pandas as pd
import tempfile, shutil
import os
import numpy as np
import pathlib
import hashlib
import pytest
from generate_test_data import make_test_data
from pciSeq.src.core.utils import get_out_dir
from pandas.testing import assert_frame_equal


# set up the logger
logger_setup()


def calculate_checksum(str_path):
    """
    Calculate SHA256 hash/checksum of a file

    Args:
        path (str): Path to file

    Returns:
        str: checksum as hexadecimal string

    From :
    https://github.com/metoppv/improver/blob/master/improver_tests/acceptance/acceptance.py
    """
    hasher = hashlib.sha256()
    path = pathlib.Path(str_path)
    with open(path, mode="rb") as kgo_file:
        while True:
            # read 1 megabyte binary chunks from file and feed them to hasher
            kgo_chunk = kgo_file.read(2 ** 20)
            if not kgo_chunk:
                break
            hasher.update(kgo_chunk)
    checksum = hasher.hexdigest()
    return checksum


def test_tmpdir():
    out_dir = os.path.join(get_out_dir(), 'tests')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


# def read_dataframe():
#     """
#     This is for the test at the bottom, needs a lot more work
#     """
#     cellData_path = os.path.join(get_out_dir(), 'pciSeq', 'data', 'cellData.tsv')
#     cellData = pd.read_csv(cellData_path, sep='\t')
#     cellData[['Cell_Num', 'X', 'Y']].round(6)
#     cellData['CellGeneCount'] = cellData['CellGeneCount'].apply(lambda x: np.round(np.array(eval(x)), 3))
#     cellData['Prob'] = cellData['Prob'].apply(lambda x: np.round(np.array(eval(x)), 3))
#     return cellData


bbox = [
    (4238, 364),  # bottomleft, [x0, y0]
    (5160, 933)  # topright, [x1, y1]
]

spots, image_label = make_test_data(test_tmpdir(), bbox=None)
coo = coo_matrix(image_label)

scRNAseq = pd.read_csv('test_scRNAseq.csv').set_index('gene_name')

# main task
# _opts = {'max_iter': 10}
opts = {'save_data': True,
         'launch_viewer': False,
         'launch_diagnostics': False,
         'output_path': [test_tmpdir()]
        }

test_list = [
    spots,
    coo,
    scRNAseq,
    opts
]

expected_list = ['51271241d3dc9d922c06d3f4e81ee57a7415da061bb80ce7b49ff97160f5d013']


@pytest.mark.parametrize('test_data, expected', [
    (test_list, expected_list)
])
def test_app(test_data, expected):
    expected_csum = expected[0]
    cellData, geneData = fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)

    csum_pickle = calculate_checksum(os.path.join(pathlib.Path(test_tmpdir()), 'pciSeq', 'data', 'debug', 'pciSeq.pickle'))
    print(csum_pickle)
    assert csum_pickle == expected_csum

    # Clean the temporary dir
    shutil.rmtree(test_tmpdir())


# THE TEST BELOW NEEDS MORE WORK
# test_cellData_df = read_dataframe()
# expected_cellData_df = read_dataframe()


# @pytest.mark.parametrize('test_cellData, expected_cellData', [
#     (test_cellData_df, expected_cellData_df)
# ])
# def test_celldata_df(test_cellData, expected_cellData):
#     ## This needs more work. it is not finished!
#     assert assert_frame_equal(test_cellData, expected_cellData)