import os
import tempfile, shutil
import pandas as pd
import numpy as np
import fastremap
import skimage.measure as skmeas
from scipy.sparse import coo_matrix, save_npz, load_npz
import logging

gtd_logger = logging.getLogger(__name__)


def make_test_data(test_data_dir, bbox=None):
    image_label = read_segmentation()
    if bbox is None:
        bbox = [
            (0, 0),                     # bottomleft, [x0, y0]
            image_label.shape[::-1]     # topright, [x1, y1]
        ]
    props = skmeas.regionprops_table(image_label, properties=['label', 'centroid'])
    props_df = pd.DataFrame(props).rename(columns={'centroid-0': 'y', 'centroid-1': 'x'})

    image_label, props_df = clip_image(image_label.copy(), props_df.copy(), bbox)
    coo_file = os.path.join(test_data_dir, 'test_label_image_coo.npz')
    save_npz(coo_file, coo_matrix(image_label))
    gtd_logger.info("saved at: %s" % coo_file)

    spots = read_spots()
    spots = clip_spots(spots.copy(), bbox)
    spots_file = os.path.join(test_data_dir, 'test_spot_data.csv')
    spots.to_csv(spots_file, index=False)
    gtd_logger.info("saved at: %s" % spots_file)

    return spots, image_label


def read_segmentation():
    coo = load_npz(os.path.join('..', 'pciSeq', 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'))
    return coo.toarray().astype(np.int32)


def clip_image(im, im_props, bbox):
    im_props = clip_dataframe(im_props.copy(), bbox)

    # Clip the label image between the bounding box (4238, 364) and (5160, 933)
    x0, y0 = bbox[0]
    x1, y1 = bbox[1]
    image_label = im[y0:y1+1, x0:x1+1]
    image_label = fastremap.mask_except(image_label, im_props.label.tolist())
    return image_label, im_props


def read_spots():
    # Load now the spots
    spots = pd.read_csv(os.path.join('..', 'pciSeq', 'data', 'mouse', 'ca1', 'iss', 'spots.csv'))
    return spots


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
    idx_x = (df.x > x0) & (df.x < x1)
    idx_y = (df.y > y0) & (df.y < y1)

    idx = idx_x & idx_y
    return df[idx]


# def get_out_dir():
#     out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq', 'tests')
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
#     return out_dir





