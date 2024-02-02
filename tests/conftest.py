import os
import pytest
import tempfile
import pandas as pd
from scipy.sparse import load_npz
import logging

gtd_logger = logging.getLogger(__name__)


bbox = [
    (4238, 364),  # bottomleft, [x0, y0]
    (5160, 933)  # topright, [x1, y1]
]


def pytest_configure():
    """
    https://docs.pytest.org/en/7.1.x/deprecations.html#pytest-namespace
    """
    pytest.fspots = None
    pytest.fcells = None
    pytest.fscData = None


def get_out_dir():
    out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq', 'tests')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


@pytest.fixture
def read_demo_data(bbox=None):
    path_str = os.path.join('..', 'pciSeq', 'data', 'mouse', 'ca1')
    spots = pd.read_csv(os.path.join(path_str, 'iss', 'spots.csv'))
    label_image = load_npz(os.path.join(path_str, 'segmentation', 'label_image.coo.npz'))

    scData = pd.read_csv('test_scRNAseq.csv').set_index('gene_name')
    if bbox is not None:
        spots, label_image = clip_data(spots.copy(), label_image.copy, bbox)
    return spots, label_image, scData


def clip_data(spots, img, bbox):
    spots_out = clip_spots(spots, bbox)
    img_out = clip_label_image(img, bbox)
    return spots_out, img_out


def clip_spots(spots, bbox):
    spots = clip_dataframe(spots.copy(), bbox)

    # Save the test spots
    x0, y0 = bbox[0]
    spots.x = spots.x - x0
    spots.y = spots.y - y0
    return spots


def clip_dataframe(df, bbox):
    x0, y0 = bbox[0]
    x1, y1 = bbox[1]
    idx_x = (df.x >= x0) & (df.x <= x1)
    idx_y = (df.y >= y0) & (df.y <= y1)

    idx = idx_x & idx_y
    return df[idx]


def clip_label_image(im, bbox):
    x0, y0 = bbox[0]
    x1, y1 = bbox[1]
    return im[y0:y1, x0:x1]



