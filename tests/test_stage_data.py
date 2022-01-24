import pandas as pd
import numpy as np
import cv2
from pciSeq.src.preprocess.spot_labels import inside_cell
import pytest

@pytest.fixture(scope='module')
def label_image():
    n = 5
    m = 25
    label_image = np.arange(n * n, dtype=np.uint16).reshape((n, n))
    label_image = cv2.resize(label_image, (m, m), interpolation=cv2.INTER_NEAREST)
    yield label_image


def test_inside_cell(label_image):
    row, col = label_image.shape
    for i in range(row):
        for j in range(col):
            point = np.array([[i,j]]).T
            v = label_image[i, j]
            assert np.all(inside_cell(label_image, point) == np.array(v))

def test_AssertTrue():
    assert True