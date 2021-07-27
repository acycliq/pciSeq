"""
Functions for get the label image given the coords of the polygons' outer rings
"""
import pandas as pd
import pciSeq.src.preprocess.vizgen.config as config
from pciSeq.src.preprocess.vizgen.utils import transformation
from scipy.sparse import coo_matrix, csr_matrix
from multiprocessing import Pool, cpu_count, freeze_support
from functools import partial
import os
from skimage.measure import regionprops
import numpy as np
from PIL import Image, ImageDraw
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


# def transformation(cgf):
#     with open(cgf['manifest']) as f:
#         settings = json.load(f)
#
#     # bounding box in microns
#     bbox = {'x0': settings['bbox_microns'][0],
#             'x1': settings['bbox_microns'][2],
#             'y0': settings['bbox_microns'][1],
#             'y1': settings['bbox_microns'][3]}
#
#     # image width and height in pixels
#     img = {'width': settings['mosaic_width_pixels'],
#            'height': settings['mosaic_height_pixels']}
#
#     # Affine transformation: a set of coefficients a, b, c, d for transforming
#     # a point of a form (x, y) into (a*x + b, c*y + d)
#     a = img['width'] / (bbox['x1'] - bbox['x0'])
#     b = -1 * img['width'] / (bbox['x1'] - bbox['x0']) * bbox['x0']
#     c = img['height'] / (bbox['y1'] - bbox['y0'])
#     d = -1 * img['height'] / (bbox['y1'] - bbox['y0']) * bbox['y0']
#
#     tx = lambda x: a * x + b
#     ty = lambda y: c * y + d
#     return tx, ty, img, bbox


# def main():
#     cfg = config.DEFAULT
#     boundaries_file = "../dashboard/data/cellBoundaries.tsv"
#     boundaries = pd.read_csv(boundaries_file, sep='\t')
#     tx, ty, img, _ = transformation(cfg)
#
#     rr = []
#     cc = []
#     data = []
#     # boundaries = boundaries.head()
#     for index, row in boundaries.iterrows():
#         poly = row['coords']
#         poly = np.array(json.loads(poly))
#         x_px = tx(poly[:,0]).astype(np.uint32)
#         y_px = ty(poly[:,1]).astype(np.uint32)
#         _rr, _cc = polygon(x_px, y_px)
#         _data = np.ones(_cc.shape) * (index+1)
#         rr.extend(_rr.astype(np.uint32))
#         cc.extend(_cc.astype(np.uint32))
#         data.extend(_data)
#         logger.info('Doing cell %d out of %d' % (index, boundaries.shape[0]))
#
#
#     # rr = np.concatenate(rr)
#     # cc = np.concatenate(cc)
#     # data = np.concatenate(data)
#
#     coo = coo_matrix((data, (rr, cc)), dtype=np.uint32)
#     # coo = coo_matrix((data, (rr, cc)), shape=(img['height'], img['width']))
#     print('Label image has shape: ',  coo.shape)
#     print('ok')


def get_mask(boundaries, label):
    """
    returns the mask from a polygon's boundaries coords
    :return:
    """
    img = Image.new('L', (2000, 2000), 0)
    ImageDraw.Draw(img).polygon(boundaries, outline=1, fill=1)
    mask = np.array(img).astype(np.int64) * label
    return mask


def worker(el, dic):
    """One worker started per CPU. Receives the label image once and a list of the labels to look for."""
    pid = os.getpid()
    # print(f'Worker pid: {pid}, processing cell id: {el[1]}')

    # cfg = config.DEFAULT
    # boundaries_file = "../dashboard/data/cellBoundaries.tsv"
    # boundaries_df = pd.read_csv(boundaries_file, sep='\t')


    _dic = {}
    cell_label = el[1]
    if cell_label % 1000 == 0:
        logger.info('Worker pid: %d processing cell id: %d' % (pid, cell_label))

    # 1. get the boundaries coords
    poly = np.array(el[2])
    # poly = np.array(json.loads(poly))

    # # 2. Convert micron to pixel coords
    # x_px = tx(poly[:,0]).astype(np.uint32)
    # y_px = ty(poly[:,1]).astype(np.uint32)

    x_px = poly[:, 0]
    y_px = poly[:, 1]

    # 3. shift the coords
    offset_x = x_px.min()
    offset_y = y_px.min()
    x = x_px - offset_x
    y = y_px - offset_y

    # 4. fill now the polygon, label it. Make a mask
    boundaries = list(zip(x, y))
    mask = get_mask(boundaries, cell_label)
    coo_mask = coo_matrix(mask)  # coo_matrix will not have duplicated coords since it is derived from the 2d array

    # 5. make now a sparse matrix of the label image. Do not forget to shift back the coords
    r = coo_mask.row + offset_y
    c = coo_mask.col + offset_x
    d = coo_mask.data

    assert np.all(np.unique(d) == cell_label)
    label_image[r, c] = d

    return dic




def get_label_image_par(list_of_lists, cfg):
    freeze_support()

    tx, ty, img, _ = transformation(cfg)
    label_image = np.zeros((img['height'], img['width']), dtype=np.uint32)

    # Get number of cores and split labels across that many workers
    processes = cpu_count()

    print(f'Using {processes} processes')

    # Chunk up the labels across the processes
    # chunks = np.array_split(list_of_lists, processes)
    chunks = np.array_split(np.array(list_of_lists, dtype="object"), processes)

    dic = {}
    # Map the labels across the processes
    pool = Pool(processes=processes)
    res = pool.map(partial(worker, label_image=label_image), list_of_lists)
    pool.close()
    pool.join()

    print('ok')
    return label_image, None


def get_label_image(boundaries_df, cfg):
    # cfg = config.DEFAULT
    # boundaries_file = "../dashboard/data/cellBoundaries.tsv"
    # boundaries_df = pd.read_csv(boundaries_file, sep='\t')
    tx, ty, img, _ = transformation(cfg)

    # label_image = csr_matrix(([], ([], [])), shape=(img['height'], img['width'])).tolil().astype(np.uint32)
    label_image = np.zeros([img['height'], img['width']]).astype(np.uint32)
    cell_props_list = []
    for index, row in boundaries_df.iterrows():
        cell_label = row['cell_label']
        if index % 1000 == 0:
            logger.info('Doing cell id: %d from a total %d' % (cell_label, boundaries_df.shape[0]))
        assert int(index + 1) == cell_label, 'Are these miss-alinged?'

        # 1. get the boundaries coords
        poly = np.array(row['cell_boundaries'])
        # poly = np.array(json.loads(poly))

        # # 2. Convert micron to pixel coords
        # x_px = tx(poly[:,0]).astype(np.uint32)
        # y_px = ty(poly[:,1]).astype(np.uint32)

        x_px = poly[:, 0]
        y_px = poly[:, 1]

        # 3. shift the coords
        offset_x = x_px.min()
        offset_y = y_px.min()
        x = x_px - offset_x
        y = y_px - offset_y

        # 4. fill now the polygon, label it. Make a mask
        boundaries = list(zip(x, y))
        mask = get_mask(boundaries, cell_label)

        # 5. Get area area and cell centroid
        props = regionprops(mask)[0]
        _cell_props = pd.DataFrame({'cell_label': [cell_label],
                                    'cell_key': row['cell_key'],
                                    'area': [props.area],
                                    'centroid_x': [props.centroid[1] + offset_x],
                                    'centroid_y': [props.centroid[0] + offset_y]})

        coo_mask = coo_matrix(mask)  # coo_matrix will not have duplicated coords since it is derived from the 2d array

        # 5. make now a sparse matrix of the label image. Do not forget to shift back the coords
        r = coo_mask.row + offset_y
        c = coo_mask.col + offset_x
        d = coo_mask.data
        # _coo = coo_matrix((mask.data, (mask.row+offset_y, mask.col+offset_x)), shape=(img['height'], img['width']))
        label_image[r, c] = d  # If further down the loop the same (r, c) appears then the label_array will
        #                        keep the most recent value

        # append now the cell props
        cell_props_list.append(_cell_props)

    # concatenate all cell_props dataframes
    cell_props = pd.concat(cell_props_list, ignore_index=True)
    coo_label_image = coo_matrix(label_image)
    return coo_label_image, cell_props


if __name__ == "__main__":
    cfg = config.DEFAULT
    boundaries_file = "../dashboard/data/cellBoundaries.tsv"
    boundaries_df = pd.read_csv(boundaries_file, sep='\t')
    label_image, cell_props = get_label_image(boundaries_df, cfg)
    logger.info(label_image.shape)
