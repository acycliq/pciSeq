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
from pciSeq.app import fit
from pciSeq.src.core.main import VarBayes
import pciSeq.config as config
from pciSeq.app import parse_args, init, validate, stage_data, cell_type
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


def tmpdir():
    out_dir = os.path.join(get_out_dir(), 'tests')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


def test_2(get_test_data):
    spots = get_test_data[0]
    coo = get_test_data[1]
    scData = get_test_data[2]

    spots, coo, scRNAseq, opts = parse_args(spots=spots, coo=coo, scRNAseq=scData, opts=config.DEFAULT)

    # 2. get the hyperparameters
    cfg = init(opts)

    # 3. validate inputs
    spots = validate(spots.copy(), coo, scRNAseq)   # Spots might get mutated here. Genes not found
                                                    # in the single cell data will be removed.

    _cells, cellBoundaries, _spots = stage_data(spots, coo)

    cellData, geneData, varBayes = cell_type(_cells, _spots, scRNAseq, cfg)


    assert len(get_test_data) == 3


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

spots, image_label = make_test_data(tmpdir(), bbox=None)
coo = coo_matrix(image_label)

scRNAseq = pd.read_csv('test_scRNAseq.csv').set_index('gene_name')

# main task
# _opts = {'max_iter': 10}
opts = {'save_data': True,
         'launch_viewer': False,
         'launch_diagnostics': False,
         'output_path': [tmpdir()]
        }

test_list = [
    spots,
    coo,
    scRNAseq,
    opts
]

expected_list = ['364e544d2c837902645083b8c2b9298fef2987a24c85160ec61c50d2e6f7ffce']


@pytest.mark.parametrize('test_data, expected', [
    (test_list, expected_list)
])
def test_app(test_data, expected):
    expected_csum = expected[0]
    cellData, geneData = fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)

    csum_pickle = calculate_checksum(os.path.join(pathlib.Path(tmpdir()), 'pciSeq', 'data', 'debug', 'pciSeq.pickle'))
    print(csum_pickle)
    assert csum_pickle == expected_csum

    # Clean the temporary dir
    shutil.rmtree(tmpdir())




# THE TEST BELOW NEEDS MORE WORK
# test_cellData_df = read_dataframe()
# expected_cellData_df = read_dataframe()


# @pytest.mark.parametrize('test_cellData, expected_cellData', [
#     (test_cellData_df, expected_cellData_df)
# ])
# def test_celldata_df(test_cellData, expected_cellData):
#     ## This needs more work. it is not finished!
#     assert assert_frame_equal(test_cellData, expected_cellData)