import os
import pandas as pd
import numpy as np
import tempfile
from typing import Tuple
from scipy.sparse import coo_matrix, save_npz, load_npz
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.viewer.utils import splitter_mb
from PIL import Image, ImageOps, ImageDraw
from pciSeq import config
import logging

logger = logging.getLogger(__name__)

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


def fit(iss_spots: pd.DataFrame, coo: coo_matrix, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for pciSeq.

    Parameters
    ----------
    iss_spots : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Gene, dtype: string, The gene name
            Name: x, dtype: int64, X-axis coordinate of the spot
            Name: y, dtype: int64, Y-axis coordinate of the spot

    coo : scipy.sparse.coo_matrix
        A label image array as a coo_matrix datatype. The label denote
        which cell the corresponding pixel 'belongs' to. If label is
        zero, the pixel is on the background

    scRNAseq : pandas.DataFrame
        Index:
            The gene name
        Columns:
            The column headers are the cell classes and the data are uint32

    opts : dictionary (Optional)
        A dictionary to pass-in user-defined hyperparameter values. They override the default
        values as these are set by the config.py file. For example to exclude genes Npy and
        Vip you can define opts as:
            opts = {'exclude_genes': ['Npy', 'Vip']}
        and pass that dict to the fit function as the last argument

    Returns
    ------
    cellData : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Cell_Num, dtype: int64, The label of the cell
            Name: X, dtype: float64, X-axis coordinate of the cell centroid
            Name: Y, dtype: float64, Y-axis coordinate of the cell centroid
            Name: Genenames, dtype: Object, array-like of the genes assinged to the cell
            Name: CellGeneCount, dtype: Object,array-like of the corresponding gene counts
            Name: ClassName, dtype: Object, array-like of the genes probable classes for the cell
            Name: Prob, dtype: Object, array-like array-like of the corresponding cell class probabilities

    geneData : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Gene, dtype: string, The gene name
            Name: Gene_id, dtype: int64, The gene id, the position of the gene if all genes are sorted.
            Name: x, dtype: int64, X-axis coordinate of the spot
            Name: y, dtype: int64, Y-axis coordinate of the spot
            Name: neighbour, dtype: int64, the label of the cell which is more likely to 'raise' the spot. If zero then the spot is a misread.
            Name: neighbour_array, dtype: Object, array-like with the labels of the 4 nearest cell. The last is always the background and has label=0
            Name: neighbour_prob, dtype: Object, array-like with the prob the corresponding cell from neighbour_array has risen the spot.
    """

    try:
        scRNAseq = kwargs['scRNAseq']
    except KeyError:
        scRNAseq = None
    except Exception as err:
        raise

    try:
        opts = kwargs['opts']
    except KeyError:
        opts = None
    except Exception as err:
        raise


    # 1. get the hyperparameters
    cfg = init(opts)

    # 2. prepare the data
    logger.info(' Preprocessing data')
    _cells, cellBoundaries, _spots = stage_data(iss_spots, coo, cfg)

    # 3. cell typing
    cellData, geneData = cell_type(_cells, _spots, scRNAseq, cfg)

    # 4. save to filesystem
    if cfg['save_data']:
        write_data(cellData, geneData, cellBoundaries)

    logger.info(' Done')
    return cellData, geneData


def cell_type(_cells, _spots, scRNAseq, ini):
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)

    logger.info(' Start cell typing')
    cellData, geneData = varBayes.run()
    return cellData, geneData


def write_data(cellData, geneData, cellBoundaries):
    out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ellipsoidBorders = cellData[['Cell_Num', 'ellipsoid_border']]
    ellipsoidBorders = ellipsoidBorders.rename(columns={'Cell_Num': 'cell_id', 'ellipsoid_border': 'coords'})

    cellData = cellData.drop(['ellipsoid_border'], axis=1)
    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))

    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))

    cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellBoundaries.tsv')))

    ellipsoidBorders.to_csv(os.path.join(out_dir, 'ellipsoidBorders.tsv'), sep='\t', index=False)
    logger.info('Saved at %s' % (os.path.join(out_dir, 'ellipsoidBorders.tsv')))

    # # Write to the disk as tsv of 99MB each
    # splitter_mb(cellData, os.path.join(out_dir, 'cellData'), 99)
    # splitter_mb(geneData, os.path.join(out_dir, 'geneData'), 99)
    # splitter_mb(cellBoundaries, os.path.join(out_dir, 'cellBoundaries'), 99)
    # splitter_mb(ellipsoidBorders, os.path.join(out_dir, 'ellipsoidBorders'), 99)


def init(opts):
    """
    Reads the opts dict and if not None, it will override the default parameter value by
    the value that the dictionary key points to.
    If opts is None, then the default values as these specified in the config.py file
    are used without any change.
    """
    cfg = config.DEFAULT
    if opts is not None:
        default_items = set(cfg.keys())
        user_items = set(opts.keys())
        assert user_items.issubset(default_items), ('Options passed-in should be a dict with keys: %s ' % default_items)
        for item in opts.items():
            if isinstance(item[1], (int, float, list)) or isinstance(item[1](1), np.floating):
                val = item[1]
            # elif isinstance(item[1], list):
            #     val = item[1]
            else:
                raise TypeError("Only integers, floats and lists are allowed")
            cfg[item[0]] = val
            logger.info(' %s is set to %s' % (item[0], cfg[item[0]]))
    return cfg


def downscale(spots, img, ppm):
    spots.x = spots.x/ppm
    spots.y = spots.y / ppm

    # z, h, w = img.shape
    # col = np.ceil(h/ppm)
    # row = np.ceil(w/ppm)
    # out = skimage.transform.resize(img, (z, col, row))

    d = int(np.floor(ppm))
    out = None
    if img is not None:
        _img = img[:, ::d, ::d]
        out = [coo_matrix(_img[i, :, :]) for i in range(len(_img))]
    return spots, out


def expand_z(label_image, spots, ppm):
    """
    expands the z-dimension so that all X,Y,Z have the same resolution.
    Also, it scales up the z-coords of the spots
    """
    z, h, w = label_image.shape

    k = int(ppm)
    depth = z * k
    temp_img = np.zeros([depth, h, w])

    coo = []
    for i in range(temp_img.shape[0]):
        j = divmod(i, k)[0]
        _img = np.array(Image.fromarray(label_image[j]).resize((width, height), Image.NEAREST), dtype=np.uint32)
        coo.append(coo_matrix(_img))

    spots.z = spots.z * ppm
    return coo, spots


def truncate_data(label_image, spots, i,  j):
    spots_out = spots.copy()
    label_image = label_image[i:j, :, :]
    spots_out = spots_out[(spots_out.z >= i) & (spots_out.z < j)]
    spots_out.z = spots_out.z.values - i
    return label_image, spots_out




if __name__ == "__main__":
    ppm = 6.0121  # pixels per micron
    width = 5865  # width of the original image
    height = 7705 # length of the original image
    z_start = 35
    z_end = 46

    # # read some demo data
    _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'B2A3', 'spots.csv'))
    # _coo = load_npz(os.path.join(ROOT_DIR, 'data', 'B2A3', 'ca1', 'segmentation', 'label_image.coo.npz'))
    label_image = np.load(os.path.join(ROOT_DIR, 'data', 'B2A3', 'B2A3_label_image.npz'))
    label_image = label_image['arr_0']  # this is already downscaled, ppm = 6.0121

    label_image, _iss_spots = truncate_data(label_image, _iss_spots, z_start,  z_end)

    # _coo = [coo_matrix(d) for d in label_image]
    _coo, _iss_spots = expand_z(label_image, _iss_spots, ppm)

    # read some demo data
    # _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'tugrul', 'TO123_S1', 'spots_shifted.csv'))
    # _coo = load_npz(os.path.join(ROOT_DIR, 'data', 'tugrul', 'TO123_S1', 'label_image.coo.npz'))

    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(float).astype(np.uint32)

    # main task
    # _opts = {'max_iter': 10}
    # my_label_image = np.stack([_coo.toarray() for i in range(10)])
    # _iss_spots, _ = downscale(_iss_spots, None, ppm)

    fit(_iss_spots, _coo, scRNAseq=_scRNAseq, opts={'save_data': True})

