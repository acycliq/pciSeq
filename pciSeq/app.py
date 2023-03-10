import os
import pandas as pd
import numpy as np
import sysconfig
import shutil
import json
import tempfile
import pickle
from typing import Tuple
from scipy.sparse import coo_matrix, save_npz, load_npz
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.cell_call.utils import get_out_dir
from pciSeq.src.preprocess.utils import get_img_shape
from pciSeq.src.viewer.run_flask import flask_app_start
from pciSeq import config
from pciSeq.src.cell_call.log_config import attach_to_log, logger

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


def fit(iss_spots: pd.DataFrame, coo: coo_matrix, scRNAseq: pd.DataFrame, opts: dict = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
            Name: Prob, dtype: Object, array-like of the corresponding cell class probabilities

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

    # 1. get the hyperparameters
    cfg = init(opts)

    # 2. validate inputs
    validate(iss_spots, scRNAseq)

    # 3. prepare the data
    logger.info(' Preprocessing data')
    _cells, cellBoundaries, _spots = stage_data(iss_spots, coo)

    # 4. cell typing
    cellData, geneData, varBayes = cell_type(_cells, _spots, scRNAseq, cfg)

    # 5. save to the filesystem
    if (cfg['save_data'] and varBayes.has_converged) or cfg['launch_viewer']:
        write_data(cellData, geneData, cellBoundaries, varBayes, cfg['output_path'])

    if cfg['launch_viewer']:
        [h, w] = get_img_shape(coo)
        dst = copy_viewer_code(cfg)
        make_config_js(dst, w, h)
        flask_app_start(dst)

    logger.info(' Done')
    return cellData, geneData


def cell_type(_cells, _spots, scRNAseq, ini):
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)

    logger.info(' Start cell typing')
    cellData, geneData = varBayes.run()
    return cellData, geneData, varBayes


def write_data(cellData, geneData, cellBoundaries, varBayes, path):
    if path[0] == 'default':
        out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq', 'data')
    else:
        out_dir = path[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))

    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))

    cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellBoundaries.tsv')))

    ## commenting this as there is no need to save the pickle file at the moment
    # with open(os.path.join(out_dir, 'pciSeq.pickle'), 'wb') as outf:
    #     pickle.dump(varBayes, outf)
    #     logger.info(' Saved at %s' % os.path.join(out_dir, 'pciSeq.pickle'))


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


def validate(spots, sc):
    assert isinstance(spots, pd.DataFrame) and set(spots.columns) == {'Gene', 'x', 'y'}, \
        "Spots should be passed-in to the fit() method as a dataframe with columns ['Gene', 'x', 'y']"

    assert isinstance(sc, pd.DataFrame), "Single cell data should be passed-in to the fit() method as a dataframe"


def make_config_base(dst):
    cellData_tsv = os.path.join(dst, 'data', 'cellData.tsv')
    geneData_tsv = os.path.join(dst, 'data', 'geneData.tsv')

    cellData_dict = {"mediaLink": "../../data/cellData.tsv", "size": str(os.path.getsize(cellData_tsv))}
    geneData_dict = {"mediaLink": "../../data/geneData.tsv", "size": str(os.path.getsize(geneData_tsv))}

    return {
        'cellData': cellData_dict,
        'geneData': geneData_dict,
    }


def make_config_js(dst, w, h):
    appDict = make_config_base(dst)
    cellBoundaries_tsv = os.path.join(dst, 'data', 'cellBoundaries.tsv')
    cellBoundaries_dict = {"mediaLink": "../../data/cellBoundaries.tsv", "size": str(os.path.getsize(cellBoundaries_tsv))}
    roi_dict = {"x0": 0, "x1": w, "y0": 0, "y1": h}
    appDict['cellBoundaries'] = cellBoundaries_dict
    appDict['roi'] = roi_dict
    appDict['zoomLevels'] = 10
    appDict['tiles'] = "https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg"

    config_str = "// NOTES: \n" \
                 "// 1. paths are with respect to the location of 'streaming-tsv-parser.js \n" \
                 "// 2. roi is the image size in pixels. Leave x0 and y0 at zero and set x1 to the width and y1 to the height \n" \
                 "// 3. tiles should point to the folder that keeps your pyramid of tiles. If you do not have that just \n" \
                 "//    change the link to a blind one (change the jpg extension for example). The viewer should work \n" \
                 "//    without the dapi background though \n" \
                 "// 4. size is the tsv size in bytes. I use os.path.getsize() to get it. Not crucial if you \n" \
                 "//    dont get it right, ie the full tsv will still be parsed despite this being wrong. It \n" \
                 "//    is used by the loading page piecharts to calc how far we are \n" \
                 "// 5. Leave zoomLevels to 10 \n" \
                 " function config() { return %s }" % json.dumps(appDict)
    config = os.path.join(dst, 'viewer', 'js', 'config.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    logger.info(' viewer config saved at %s' % config)


def copy_viewer_code(cfg):
    site_packages_dir = sysconfig.get_path('purelib')
    pciSeq_dir = os.path.join(site_packages_dir, 'pciSeq')
    dim = '2D'
    src = os.path.join(pciSeq_dir, 'static', dim)
    dst = get_out_dir(cfg['output_path'], '')

    shutil.copytree(src, dst, dirs_exist_ok=True)
    logger.info(' viewer code (%s) copied from %s to %s' % (dim, src, dst))
    return dst


if __name__ == "__main__":
    # set up the logger
    attach_to_log()

    # read some demo data
    _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'iss', 'spots.csv'))
    _coo = load_npz(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'))

    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(float).astype(np.uint32)

    # main task
    # _opts = {'max_iter': 10}
    fit(_iss_spots, _coo, _scRNAseq, {'save_data': True})

