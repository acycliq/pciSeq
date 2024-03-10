import os
import pytest
import pathlib
import hashlib
import numpy as np
import pciSeq.config as config
from pciSeq.src.core.main import VarBayes
from pciSeq.src.core.utils import get_out_dir
from pciSeq.app import parse_args, validate, stage_data, init, fit
from pciSeq.src.diagnostics.launch_diagnostics import launch_dashboard
import logging


def calculate_checksum(str_path):
    """
    Calculate SHA256 hash/checksum of a file
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
    scData = read_demo_data[2]

    with pytest.raises(AssertionError) as excinfo:
        validate(coo, coo, scData, config.DEFAULT)
    assert str(excinfo.value) == ("Spots should be passed-in to the fit() method as "
                                  "a dataframe with columns ['Gene', 'x', 'y']")


# -------------------------------------------------------------------
# first element is the expected area.sum
# second element is the expected max label
# third element is the expected column-wize sum of the spots dataframe
expected_list = [
    3721912.0,
    3481,
    [231498829.0,150952363.0,75868924.0,151411232.0,82555688.0]
    # [231498829.00000, 150952363.00000, 75868924.00000, 151411235.74292, 82555694.21166]
]

231498829.0,150952363.0,75868924.0,151411232.0,82555688.0


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
expected_iter_delta = [1.0, 0.8998741, 0.67078424, 0.47034895,
                       0.6274648, 0.5809673, 0.60502154, 0.28165847,
                       0.26288757, 0.19117153, 0.5582824, 0.1904957,
                       0.24826217, 0.28257173, 0.068271875, 0.1609922,
                       0.58281434, 0.16718799, 0.14889961, 0.1229637,
                       0.089731336, 0.044867665, 0.06780699, 0.10764623,
                       0.19164285, 0.16254145, 0.061941266, 0.044866115,
                       0.031517506, 0.019638658]


@pytest.mark.parametrize('filename, expected', [
    ('read_demo_data', expected_iter_delta)
])
def test_varBayes(filename, expected, request):
    read_demo_data = request.getfixturevalue(filename)
    spots = read_demo_data[0]
    coo = read_demo_data[1]
    scData = read_demo_data[2]

    cfg = init({'launch_viewer': True,
                'launch_diagnostics': True,
                'max_iter': 30,
                })
    validate(spots, coo, scData, cfg)

    if cfg['launch_diagnostics'] and cfg['is_redis_running']:
        launch_dashboard()

    # prepare the data
    _cells, cellBoundaries, _spots = stage_data(spots, coo)

    # cell typing
    varBayes = VarBayes(_cells, _spots, scData, cfg)
    cellData, geneData = varBayes.run()

    arr_1 = np.array(varBayes.iter_delta, dtype=np.float32).round(11)
    arr_2 = np.array(expected_iter_delta[: cfg['max_iter']], dtype=np.float32).round(11)
    assert np.all(arr_1 == arr_2)


