import pandas as pd
import numpy as np
import cv2
from distutils import dir_util
import os
from scipy.sparse import load_npz
from pciSeq.src.preprocess.spot_labels import inside_cell, extract_borders_dip, stage_data
from pciSeq.src.cell_call import utils
import pytest

@pytest.fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir

@pytest.fixture
def ground_truth_spots(datadir):
    return pd.read_csv(datadir.join('results_spots.csv'))

@pytest.fixture(scope='module')
def dummy_label_image():
    n = 5
    m = 25
    label_image = np.arange(n * n, dtype=np.uint16).reshape((n, n))
    label_image = cv2.resize(label_image, (m, m), interpolation=cv2.INTER_NEAREST)
    yield label_image

@pytest.fixture(scope='module')
def coo_label_image():
    coo_file = utils.load_from_url(
        'https://github.com/acycliq/pciSeq/blob/dev/pciSeq/data/mouse/ca1/segmentation/label_image.coo.npz?raw=true')
    coo = load_npz(coo_file)
    yield coo

@pytest.fixture(scope='module')
def spots():
    spots_file = utils.load_from_url('https://github.com/acycliq/pciSeq/blob/dev/pciSeq/data/mouse/ca1/iss/spots.csv?raw=true')
    yield pd.read_csv(spots_file)


def test_inside_cell(spots, coo_label_image):
    inc = inside_cell(coo_label_image.tocsr(), spots)

    # background spots
    assert sum(inc == 0) == 27983

    # spots within a cell
    assert sum(inc > 0) == 44353

    # total number of spots
    assert sum(inc == 0) + sum(inc > 0) == spots.shape[0]


def test_extract_borders_dip(coo_label_image):
    cell_boundaries = extract_borders_dip(coo_label_image.toarray().astype(np.uint32), 0, 0, [0])
    assert isinstance(cell_boundaries, pd.DataFrame)

    assert cell_boundaries.label.max() == 3481

    out = 0
    for index, row in cell_boundaries.iterrows():
        out += np.array(row['coords']).sum()
    assert out == 505840577.0

def test_stage_data(spots, coo_label_image, datadir):
    cells, cell_boundaries, dots = stage_data(spots, coo_label_image)
    assert isinstance(cells, pd.DataFrame)
    assert isinstance(cell_boundaries, pd.DataFrame)
    assert isinstance(dots, pd.DataFrame)

    assert cells.label.max() == 3481.0
    assert cells.area.sum() == 3721912.0
    assert cells[['x', 'y']].values.sum() == 18715156.370488837

    out = 0
    for index, row in cell_boundaries.iterrows():
        out += np.array(row['coords']).sum()
    assert out == 505840577.0

    spots_good = pd.read_csv(datadir.join('results_spots.csv'))
    assert np.nanmax(abs(dots[['x_global', 'y_global', 'label', 'x_cell', 'y_cell']] - spots_good[['x_global', 'y_global', 'label', 'x_cell', 'y_cell']])) < 1e-10












