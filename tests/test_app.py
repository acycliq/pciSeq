from pciSeq.src.core.logger import logger_setup
from pciSeq import fit
from scipy.sparse import load_npz, coo_matrix
import pandas as pd
import tempfile, shutil
import os
import pathlib
import hashlib
import pytest
from generate_test_data import make_test_data


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


def get_out_dir():
    out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq', 'tests')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


bbox = [
    (4238, 364),  # bottomleft, [x0, y0]
    (5160, 933)  # topright, [x1, y1]
]

spots, image_label = make_test_data(get_out_dir(), bbox=None)
coo = coo_matrix(image_label)

scRNAseq = pd.read_csv('test_scRNAseq.csv', header=None, dtype=object)
scRNAseq.columns = scRNAseq.iloc[0]
scRNAseq = scRNAseq.iloc[1:, :]
scRNAseq = scRNAseq.set_index('gene_name')
scRNAseq = scRNAseq.applymap(lambda x: eval(x))

# main task
# _opts = {'max_iter': 10}
opts = {'save_data': True,
         'launch_viewer': False,
         'launch_diagnostics': False,
         'output_path': [get_out_dir()]
        }

test_list = [
    spots,
    coo,
    scRNAseq,
    opts
]

expected_list = ['2e1264e90c12aebf51717f8faab27aa95eb620175d013aa36163fff414caed08']


@pytest.mark.parametrize('test_data, expected', [
    (test_list, expected_list)
])
def test_app(test_data, expected):
    expected_csum = expected[0]
    cellData, geneData = fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)

    tmp_dir = os.path.join(get_out_dir(), 'pciSeq', 'data')
    csum_pickle = calculate_checksum(os.path.join(pathlib.Path(tmp_dir), 'debug', 'pciSeq.pickle'))
    print(csum_pickle)
    assert csum_pickle == expected_csum

    # Clean the temporary dir
    shutil.rmtree(get_out_dir())

