import os
import numpy as np
import pathlib
import hashlib
import pytest
from pciSeq.src.core.utils import get_out_dir
from pciSeq.src.core.main import VarBayes
import pciSeq.config as config
from pciSeq.app import parse_args, validate, stage_data
import logging


# set up the logger
# logger_setup()


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


# -------------------------------------------------------------------
def test_parse_args(read_demo_data):
    logging.getLogger().info('test_2')
    spots = read_demo_data[0]
    coo = read_demo_data[1]

    with pytest.raises(AssertionError) as excinfo:
        parse_args(spots)
    assert str(excinfo.value) == ('Need to provide the spots and the coo matrix as the first '
                                  'and second args to the fit() method.')

    _, _, scData, opts = parse_args(spots, coo)
    assert scData is None
    assert opts is None


# -------------------------------------------------------------------
def test_validate(read_demo_data):
    logging.getLogger().info('test_2')
    spots = read_demo_data[0]
    coo = read_demo_data[1]

    with pytest.raises(AssertionError) as excinfo:
        validate(coo, coo, config.DEFAULT)
    assert str(excinfo.value) == ("Spots should be passed-in to the fit() method as "
                                  "a dataframe with columns ['Gene', 'x', 'y']")


# -------------------------------------------------------------------
# first element is the expected area.sum
# second element is the expected max label
# third element is the expected column-wize sum of the spots dataframe
expected_list = [
    3721912.0,
    3481,
    [231498829.00000, 150952363.00000, 75868924.00000, 151411235.74292, 82555694.21166]
]


@pytest.mark.parametrize('filename, expected', [
    ('read_demo_data', expected_list)
])
def test_stage_data(filename, expected, request):
    read_demo_data = request.getfixturevalue(filename)
    spots = read_demo_data[0]
    coo = read_demo_data[1]
    cells, cell_boundaries, spots = stage_data(spots, coo)
    pytest.fspots = spots
    pytest.fcells = cells
    assert len(cells.label) == len(np.unique(cells.label))
    assert cells.area.sum() == expected_list[0]
    assert cells.label.max() == expected_list[1]
    assert np.all(
        spots[['x_global', 'y_global', 'label', 'x_cell', 'y_cell']].sum().round(5).values == expected_list[2])


# -------------------------------------------------------------------
#  Expected max diff of the spot-to-cell probability between two successive iterations
expected_iter_delta = [
    1.0, 0.8998741535450602, 0.6707994416612973, 0.4703557918803584, 0.6274678470436773, 0.5809977507020871,
    0.6050537780421158, 0.281664156778341, 0.26291180844419526, 0.1911332823201044, 0.5582350505154667,
    0.1904928001599277, 0.2482987581795233, 0.2825277025986648, 0.06827397193459617, 0.16097980010421764,
    0.5828276333586281, 0.16734418317139443, 0.14889440752498986, 0.1229702319182493, 0.08985526872588012,
    0.044863509416499414, 0.0678009601806856, 0.10751329267974208, 0.19152470052524656, 0.16275473631515358,
    0.06207517767135862, 0.04486396122645342, 0.0315151843141393, 0.01963733267979928
]


@pytest.mark.parametrize('filename, expected', [
    ('read_demo_data', expected_iter_delta)
])
def test_varBayes(filename, expected, request):
    read_demo_data = request.getfixturevalue(filename)
    scData = read_demo_data[2]
    opts = config.DEFAULT
    opts['launch_viewer'] = True
    varBayes = VarBayes(pytest.fcells, pytest.fspots, scData, opts)
    varBayes.run()

    arr_1 = np.array(varBayes.iter_delta).round(11)
    arr_2 = np.array(expected_iter_delta).round(11)
    assert np.all(arr_1 == arr_2)


