import os
import glob
import csv
import numpy as np
import pandas as pd
import math
from multiprocessing import Pool, cpu_count, freeze_support
from functools import partial
import pyvips
import numba
import json
# import dropbox
import logging
import pciSeq.src.preprocess.vizgen.config as config

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def _get_file(OUT_DIR, n, header_line):
    filename = os.path.basename(OUT_DIR).split('.')[0]
    file = os.path.join(OUT_DIR, filename + '_%d.%s' % (n, 'tsv'))
    handle = open(file, "a", newline='', encoding='utf-8')
    write = csv.writer(handle, delimiter='\t')
    write.writerow(header_line)
    return file, handle


def splitter_mb(df, dir_path, mb_size):
    """ Splits a text file in (almost) equally sized parts on the disk. Assumes that there is a header in the first line
    :param filepath: The path of the text file to be broken up into smaller files
    :param mb_size: size in MB of each chunk
    :return:
    """
    # OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    else:
        files = glob.glob(dir_path + '/*.*')
        for f in files:
            os.remove(f)

    n = 0
    header_line = df.columns.tolist()
    # header_line = next(handle)[1].tolist()
    file_out, handle_out = _get_file(dir_path, n, header_line)
    # data_row = next(handle)[1].tolist()
    for index, row in df.iterrows():
        row = row.tolist()
        size = os.stat(file_out).st_size
        if size > mb_size*1024*1024:
            logger.info('saved %s with file size %4.3f MB' % (file_out, size/(1024*1024)))
            n += 1
            handle_out.close()
            file_out, handle_out = _get_file(dir_path, n, header_line)
        write = csv.writer(handle_out, delimiter='\t')
        write.writerow(row)

    # print(str(file_out) + " file size = \t" + str(size))
    logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
    handle_out.close()


def save_df(df, dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    else:
        files = glob.glob(dir_path + '/*.*')
        for f in files:
            os.remove(f)

    # Get number of cores and split labels across that many workers
    processes = cpu_count()
    logger.info(f'Using {processes} processes')

    # Chunk up the labels across the processes
    chunks = np.array_split(df, processes)

    # Map the labels across the processes
    with Pool(processes=processes) as pool:
        result = pool.map(partial(worker, OUT_DIR=dir_path), enumerate(chunks))
        pool.close()
        pool.join()


def save_df_simple(df, dir_path):
    """
    Simple convenience function to save a dataframe as a tsv.
    It creates the containing folder if it doesnt exist
    :param df:
    :param dir_path:
    :return:
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    else:
        files = glob.glob(dir_path + '/*.*')
        for f in files:
            os.remove(f)
    filename = os.path.basename(dir_path).split('.')[0]
    file_name = os.path.join(dir_path, filename + '.tsv')
    df.to_csv(file_name, sep='\t', index=False)


def worker(arg_in, OUT_DIR):
    n = arg_in[0]
    df = arg_in[1]
    filename = os.path.basename(OUT_DIR).split('.')[0]
    file = os.path.join(OUT_DIR, filename + '_%d.%s' % (n, 'tsv'))
    df.to_csv(file, index=False, sep='\t')
    # logger.info('saved %s ' % file)


def dapi_dims(cfg):
    dapi = pyvips.Image.new_from_file(cfg['dapi_tif'], access='sequential')
    img = {'width': dapi.width,
           'height': dapi.height}
    return img


def dapi_dims_from_manifest(cfg):
    with open(cfg['manifest']) as f:
        settings = json.load(f)

    # image width and height in pixels
    img = {'width': settings['mosaic_width_pixels'],
           'height': settings['mosaic_height_pixels']}

    return img


def transformation(cfg):
    # Micron to pixel transformation

    # settings = dropbox_streamer(cgf['manifest'])

    um_to_px = np.genfromtxt(cfg['micron_to_mosaic_pixel_transform'], delimiter=' ')
    assert um_to_px.shape == (3,3), 'The file %s should contain a space delimited 3-by-3 array' % cfg['micron_to_mosaic_pixel_transform']
    a = um_to_px[0, 0]
    b = um_to_px[0, 2]
    c = um_to_px[1, 1]
    d = um_to_px[1, 2]

    # img = dapi_dims(cfg)
    img = dapi_dims_from_manifest(cfg)

    # bounding box in microns
    bbox = {}
    bbox['x0'] = -1 * b/a
    bbox['x1'] = img['width']/a + bbox['x0']
    bbox['y0'] = -1 * d/c
    bbox['y1'] = img['height']/a + bbox['y0']

    # with open(cfg['manifest']) as f:
    #     settings = json.load(f)
    #
    # # bounding box in microns
    # bbox = {'x0': settings['bbox_microns'][0],
    #         'x1': settings['bbox_microns'][2],
    #         'y0': settings['bbox_microns'][1],
    #         'y1': settings['bbox_microns'][3]}
    #
    # # image width and height in pixels
    # img = {'width': settings['mosaic_width_pixels'],
    #        'height': settings['mosaic_height_pixels']}
    #
    # # Affine transformation: a set of coefficients a, b, c, d for transforming
    # # a point of a form (x, y) into (a*x + b, c*y + d)
    # a0 = img['width'] / (bbox['x1'] - bbox['x0'])
    # b0 = -1 * img['width'] / (bbox['x1'] - bbox['x0']) * bbox['x0']
    # c0 = img['height'] / (bbox['y1'] - bbox['y0'])
    # d0 = -1 * img['height'] / (bbox['y1'] - bbox['y0']) * bbox['y0']

    tx = lambda x: a * x + b
    ty = lambda y: c * y + d
    return tx, ty, img, bbox


# def dropbox_streamer(yourpath):
#     _, file_extension = os.path.splitext('/path/to/somefile.ext')
#     token = credentials.DROPBOX_TOKEN
#     dbx = dropbox.Dropbox(token)
#
#     # Relevant streamer
#     def stream_dropbox_file(path):
#         _, res = dbx.files_download(path)
#         with closing(res) as result:
#             byte_data = result.content
#             return io.BytesIO(byte_data)
#
#     # Usage
#     file_stream = stream_dropbox_file(yourpath)
#     if file_extension == '/csv':
#         out = pd.read_csv(file_stream)
#     elif file_extension == '/tsv':
#         out = pd.read_csv(file_stream, sep='\t')
#     else:
#         out = json.load(file_stream)
#     return out


def get_bbox(w, h):
    #     w = im.shape[1]
    #     h = im.shape[0]
    bottom_left = [0, 0]
    bottom_right = [w, 0]
    top_left = [0, h]
    top_right = [w, h]
    centre = [w / 2, h / 2]

    bbox = np.array([bottom_left, bottom_right, top_left, top_right, centre])
    return bbox


# 3. rotate the bounding box
def bbox_rot(w, h, R):
    offset_x = 0
    offset_y = 0
    bbox = get_bbox(w, h)
    _rot = np.array([R.dot(d).tolist() for d in bbox])

    x = _rot[:, 0].min()
    y = _rot[:, 1].min()
    if np.any(x < 0):
        offset_x = -1 * x.min()

    if np.any(y < 0):
        offset_y = -1 * y.min()
    return offset_x, offset_y


def rotate_data(points, cfg):
    """
    rotates the spots by theta_degrees
    :param my_points:
    :param theta_deg:
    :param cfg:
    :return:
    """
    theta_deg = cfg['rotation'][0]
    theta_rad = math.radians(theta_deg)

    # img = dapi_dims(cgf)
    img = dapi_dims_from_manifest(cfg)
    w = img['width']
    h = img['height']

    # points = data[['x', 'y']].values

    # 1. rotate the image
    # w, h = rotate_image(img_in, 'girl_rot.jpg', theta_deg)

    # 2. Rotation matrix to rotate the datapoints
    R = np.array([[np.cos(theta_rad), -np.sin(theta_rad)], [np.sin(theta_rad), np.cos(theta_rad)]])

    # 2 get the offsets
    offset_x, offset_y = bbox_rot(w, h, R)

    # 3 rotate the datapoints
    rot = np.dot(points, R.T)
    # logger.info('rotating finished')
    rot = rot + np.array([offset_x, offset_y])
    return rot


@numba.jit(nopython=True)
def _is_inside_sm(polygon, point):
    # From https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
    # and
    # https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1


@numba.njit(parallel=True)
def _is_inside_sm_parallel(points, polygon):
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        D[i] = _is_inside_sm(polygon,points[i])
    return D


def is_inside(points, polygon):
    # make sure polygon is closed
    if tuple(polygon[0]) != tuple(polygon[-1]):
        polygon = np.append(polygon, polygon[0, None], axis=0)
    return _is_inside_sm_parallel(points, polygon)



if __name__=="__main__":
    cfg = config.DEFAULT
    tx, ty = transformation(cfg)

    spots_path = cfg['detected_transcripts']
    chunks = pd.read_csv(spots_path, chunksize=100000)
    data = pd.concat(chunks)
    data.shape


    data_z3 = data[data.global_z == 3]
    [unique_gene_names, gene_id, counts_per_gene] = np.unique(data_z3.gene.values, return_inverse=True,
                                                              return_counts=True)
    data_z3 = data_z3[['gene', 'global_x', 'global_y']].rename(
        columns={'gene': "Gene", 'global_x': "x", 'global_y': "y"})

    data_z3['Gene_id'] = gene_id
    data_z3['neighbour'] = np.zeros(len(gene_id))
    data_z3['neighbour_array'] = [[0] for i in range(len(gene_id))]
    data_z3['neighbour_prob'] = [[1.0] for i in range(len(gene_id))]

    splitter_mb(data_z3, 'geneData', 99)
    logger.info('Done!')





