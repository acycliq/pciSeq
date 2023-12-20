import os
import pandas as pd
import numpy as np
import tempfile
import pickle
import redis
from typing import Tuple
from scipy.sparse import coo_matrix, save_npz, load_npz
from pciSeq.src.core.main import VarBayes
from pciSeq.src.core.utils import get_out_dir
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.preprocess.utils import get_img_shape
from pciSeq.src.diagnostics.utils import redis_db
from pciSeq.src.diagnostics.launch_diagnostics import launch_dashboard
from pciSeq.src.viewer.run_flask import flask_app_start
from pciSeq.src.viewer.utils import copy_viewer_code, make_config_js, make_classConfig_js
from pciSeq import config
from pciSeq.src.core.log_config import attach_to_log, logger

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


def fit(*args, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for pciSeq.

    Parameters
    ----------
    **spots : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Gene, dtype: string, The gene name
            Name: x, dtype: int64, X-axis coordinate of the spot
            Name: y, dtype: int64, Y-axis coordinate of the spot

    **coo : scipy.sparse.coo_matrix
        A label image array as a coo_matrix datatype. The label denote
        which cell the corresponding pixel 'belongs' to. If label is
        zero, the pixel is on the background

    **scRNAseq : pandas.DataFrame (Optional)
        Index:
            The gene name
        Columns:
            The column headers are the cell classes and the data are uint32

    **opts : dictionary (Optional)
        A dictionary to pass-in user-defined hyperparameter values. They override the default
        values as these are set by the config.py file. For example to exclude genes Npy and
        Vip you can define opts as:
            opts = {'exclude_genes': ['Npy', 'Vip']}
        and pass that dict to the fit function as the last argument

    *args: If 'spots' and 'coo' are not passed-in as keyword arguments then they should be provided
    as first and second positional arguments

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

    # 1. parse/check the arguments
    spots, coo, scRNAseq, opts = parse_args(*args, **kwargs)

    # 2. get the hyperparameters
    cfg = init(opts)

    # 3. validate inputs
    spots = validate(spots, coo, scRNAseq)  # Spots might be mutated here by removing entries not found
                                            # in the single cell data

    # 4. launch the diagnostics
    if cfg['launch_diagnostics'] and cfg['is_redis_running']:
        logger.info('Launching the diagnostics dashboard')
        launch_dashboard()

    # 5. prepare the data
    logger.info(' Preprocessing data')
    _cells, cellBoundaries, _spots = stage_data(spots, coo)

    # 6. cell typing
    cellData, geneData, varBayes = cell_type(_cells, _spots, scRNAseq, cfg)

    # 7. save to the filesystem
    if (cfg['save_data'] and varBayes.has_converged) or cfg['launch_viewer']:
        write_data(cellData, geneData, cellBoundaries, varBayes, cfg)

    # 8. do the viewer if needed
    if cfg['launch_viewer']:
        dst = pre_launch(cellData, coo, scRNAseq, cfg)
        flask_app_start(dst)

    logger.info(' Done')
    return cellData, geneData


def cell_type(_cells, _spots, scRNAseq, ini):
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)

    logger.info(' Start cell typing')
    cellData, geneData = varBayes.run()
    return cellData, geneData, varBayes


def write_data(cellData, geneData, cellBoundaries, varBayes, cfg):
    dst = get_out_dir(cfg['output_path'])
    out_dir = os.path.join(dst, 'data')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))

    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))

    cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellBoundaries.tsv')))

    serialise(varBayes, os.path.join(out_dir, 'debug'))


def serialise(varBayes, debug_dir):
    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)
    pickle_dst = os.path.join(debug_dir, 'pciSeq.pickle')
    with open(pickle_dst, 'wb') as outf:
        pickle.dump(varBayes, outf)
        logger.info(' Saved at %s' % pickle_dst)


def export_db_tables(out_dir, con):
    tables = con.get_db_tables()
    for table in tables:
        export_db_table(table, out_dir, con)


def export_db_table(table_name, out_dir, con):
    df = con.from_redis(table_name)
    fname = os.path.join(out_dir, table_name + '.csv')
    df.to_csv(fname, index=False)
    logger.info(' Saved at %s' % fname)


def init(opts):
    """
    Reads the opts dict and if not None, it will override the default parameter value by
    the value that the dictionary key points to.
    If opts is None, then the default values as these specified in the config.py file
    are used without any change.
    """
    cfg = config.DEFAULT
    cfg['is_redis_running'] = check_redis_server()
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


def validate(spots, coo, sc):
    assert isinstance(spots, pd.DataFrame) and set(spots.columns) == {'Gene', 'x', 'y'}, \
        "Spots should be passed-in to the fit() method as a dataframe with columns ['Gene', 'x', 'y']"

    assert isinstance(coo, coo_matrix), 'The segmentation masks should be passed-in as a coo_matrix'

    if sc is not None:
        assert isinstance(sc, pd.DataFrame), "Single cell data should be passed-in to the fit() method as a dataframe"

        if not set(spots.Gene).issubset(sc.index):
            # remove genes that cannot been found in the single cell data
            spots = purge_spots(spots, sc)

    return spots


def purge_spots(spots, sc):
    drop_spots = list(set(spots.Gene) - set(sc.index))
    logger.warning('Found %d genes that are not included in the single cell data' % len(drop_spots))
    idx = ([i for i, v in enumerate(spots.Gene) if v not in drop_spots]) # can you speed this up?? What happens if the spot df has millions of rows?
    spots = spots.iloc[idx]
    logger.warning('Removed from spots: %s' % drop_spots)
    return spots


def parse_args(*args, **kwargs):
    '''
    Do soma basic checking of the args
    '''

    spots = None
    coo = None
    scRNAseq = None
    opts = None

    if not {'spots', 'coo'}.issubset(set(kwargs)):
        try:
            assert len(args) == 2, 'Need to provide the spots and the coo matrix as the first ' \
                                   'and second args to the fit() method '
            kwargs['spots'] = args[0]
            kwargs['coo'] = args[1]
        except Exception as err:
            raise

    try:
        spots = kwargs['spots']
        coo = kwargs['coo']
    except Exception as err:
        raise

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

    return spots, coo, scRNAseq, opts


def pre_launch(cellData, coo, scRNAseq, cfg):
    '''
    Does some housekeeping, pre-flight control checking before the
    viewer is triggered

    Returns the destination folder that keeps the viewer code and
    will be served to launch the website
    '''
    [h, w] = get_img_shape(coo)
    dst = get_out_dir(cfg['output_path'])
    copy_viewer_code(cfg, dst)
    make_config_js(dst, w, h)
    if scRNAseq is None:
        label_list = [d[:] for d in cellData.ClassName.values]
        labels = [item for sublist in label_list for item in sublist]
        labels = sorted(set(labels))

        make_classConfig_js(labels, dst)
    return dst


def confirm_prompt(question):
    reply = None
    while reply not in ("", "y", "n"):
        reply = input(f"{question} (y/n): ").lower()
    return (reply in ("", "y"))


def check_redis_server():
    logger.info("check_redis_server")
    try:
        redis_db()
        return True
    except (redis.exceptions.ConnectionError, ConnectionRefusedError, OSError):
        logger.info("Redis ping failed!. Diagnostics will not be called unless redis is installed and the service is running")
        return False
        # logger.info("Redis ping failed!. Trying to install redis server")
        # if confirm_prompt("Do you want to install redis server?"):
        #     logger.info("...installing...")


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
    _opts = {'save_data': True,'launch_viewer': True, 'launch_diagnostics': True}
    fit(spots=_iss_spots, coo=_coo, scRNAseq=_scRNAseq, opts=_opts)


