import os
import pandas as pd
import numpy as np
import tempfile
from typing import Tuple
from scipy.sparse import coo_matrix, save_npz, load_npz
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.viewer.utils import splitter_mb
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
    logger.info(' Pixels per micron is set to: %f' % cfg['ppm'])

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


# def _reorder_labels(label_image):
#     """
#     rearranges the labels so that they are a sequence of integers
#     """
#     # label_image = np.stack([d.toarray().astype(np.uint16) for d in coo])
#     _, idx = np.unique(label_image.flatten(), return_inverse=True)
#     my_masks = idx.reshape(label_image.shape)
#     out_1 = [coo_matrix(d) for d in my_masks]
#     out_2 = [coo_matrix(d.astype(np.uint16)) for d in my_masks]
#     return [coo_matrix(d) for d in my_masks]

if __name__ == "__main__":

    # # read some demo data
    # _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'B2A3', 'truncated_data', 'B2A3_spots_truncated.csv'))
    # _coo = np.load(os.path.join(ROOT_DIR, 'data', 'B2A3', 'truncated_data', 'B2A3_label_image_truncated.npz'),  allow_pickle=True)['arr_0']

    # _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', '220308', 'spots_min.csv'))
    # _iss_spots = _iss_spots.assign(z=_iss_spots.z_stack * config.DEFAULT['ppm'])
    # _coo = np.load(os.path.join(ROOT_DIR, 'data',  '220308', 'coo_list.npz'),  allow_pickle=True)['arr_0']

    _iss_spots = pd.read_csv(r"E:\data\Anne\220308 50umCF seq atto425 DY520XL MS002\spots_yxz.csv")
    _iss_spots = _iss_spots.assign(z_stack=_iss_spots.z)
    _iss_spots = _iss_spots[['y', 'x', 'z_stack', 'Gene']]
    _iss_spots = _iss_spots.assign(z=_iss_spots.z_stack * config.DEFAULT['ppm'])
    _coo = np.load(r"E:\data\Anne\220308 50umCF seq atto425 DY520XL MS002\masks_2D_stiched_fullsize.npz",  allow_pickle=True)['arr_0']
    _coo = [coo_matrix(d) for d in _coo]



    # _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'B2A3', 'small_data_3d', 'small_spots.csv'))
    # logger.info('Keep Id2 spots only!')
    # _iss_spots = _iss_spots[_iss_spots.Gene == 'Id2']
    # _coo = np.load(os.path.join(ROOT_DIR, 'data', 'B2A3', 'small_data_3d', 'small_coo_from_page40.npz'), allow_pickle=True)['arr_0']

    # read some demo data

    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(float).astype(np.uint32)

    # main task
    # _opts = {'max_iter': 10}
    fit(_iss_spots, _coo, scRNAseq=_scRNAseq, opts={'save_data': True, 'z_stack_min': 18, 'z_stack_max': 43})
