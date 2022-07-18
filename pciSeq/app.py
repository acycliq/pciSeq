import os
import pandas as pd
import numpy as np
import tempfile
import pickle
from typing import Tuple
from scipy.sparse import coo_matrix, save_npz, load_npz
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq import config
from pciSeq.src.cell_call.log_config import attach_to_log, logger

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

    # 2. validate inputs
    iss_spots, coo, cfg = validate(iss_spots, coo, scRNAseq, cfg)

    # 3. prepare the data
    logger.info(' Preprocessing data')
    _cells, cellBoundaries, _spots = stage_data(iss_spots, coo, cfg)

    # 4. cell typing
    cellData, geneData, varBayes = cell_type(_cells, _spots, scRNAseq, cfg)

    # 5. save to the filesystem
    if cfg['save_data'] and varBayes.has_converged:
        write_data(cellData, geneData, cellBoundaries, varBayes, path=cfg['output_path'])

    if cfg['launch_viewer']:
        if cfg['relax_segmentation'] or cfg['is_3D']:
            pass
        else:
            pass

    varBayes.conn.close()
    logger.info(' Done')
    return cellData, geneData


def cell_type(_cells, _spots, scRNAseq, ini):
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)

    logger.info(' Start cell typing')
    cellData, geneData = varBayes.run()
    return cellData, geneData, varBayes


def write_data(cellData, geneData, cellBoundaries, varBayes, path):

    if path[0] == 'default':
        out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq')
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

    if 'ellipsoid_border' in cellData.columns:
        ellipsoidBorders = cellData[['Cell_Num', 'ellipsoid_border']]
        ellipsoidBorders = ellipsoidBorders.rename(columns={'Cell_Num': 'cell_id', 'ellipsoid_border': 'coords'})

        ellipsoidBorders.to_csv(os.path.join(out_dir, 'ellipsoidBorders.tsv'), sep='\t', index=False)
        logger.info('Saved at %s' % (os.path.join(out_dir, 'ellipsoidBorders.tsv')))

    export_db_tables(out_dir, varBayes.conn)
    # with open(os.path.join(out_dir, 'pciSeq.pickle'), 'wb') as outf:
    #     pickle.dump(varBayes, outf)
    #     logger.info(' Saved at %s' % os.path.join(out_dir, 'pciSeq.pickle'))



    # Write to the disk as tsv of 99MB each
    # splitter_mb(cellData, os.path.join(out_dir, 'cellData'), 99)
    # splitter_mb(geneData, os.path.join(out_dir, 'geneData'), 99)
    # splitter_mb(cellBoundaries, os.path.join(out_dir, 'cellBoundaries'), 99)


def export_db_tables(out_dir, con):
    str = "SELECT name FROM sqlite_schema WHERE type = 'table' ORDER BY name;"
    tables = con.execute(str).fetchall()
    for table in tables:
        export_db_table(table[0], out_dir, con)


def export_db_table(table_name, out_dir, con):
    if table_name == 'spots':
        query_str = "SELECT * FROM %s " % table_name
    else:
        query_str = "SELECT * FROM %s where has_converged = 1 " % table_name
    df = pd.read_sql_query(query_str, con)
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


def validate(spots, coo, sc, cfg):

    # check if z_stack column is missing
    if 'z_stack' not in spots.columns or not cfg['is_3D'] or not cfg['relax_segmentation']:
        spots = spots.assign(z_stack=np.zeros(spots.shape[0]))

    assert isinstance(spots, pd.DataFrame) and set(spots.columns) == {'Gene', 'x', 'y', 'z_stack'}, \
        "Spots should be passed-in to the fit() method as a dataframe with columnns ['Gene', 'x', 'y']"

    if isinstance(coo, coo_matrix):
        coo = [coo]
    assert np.all([isinstance(d, coo_matrix) for d in coo]) and isinstance(coo, list), \
        "The label image should be passed in as a coo_matrix (if you run pciSed in 2D) or as a list of coo_matrices"

    if sc is not None:
        assert isinstance(sc, pd.DataFrame), "Single cell data should be passed-in to the fit() method as a dataframe"

    if '3D:anisotropy' not in cfg or not cfg['is_3D']:
        cfg['3D:anisotropy'] = 1.0

    if '3D:from_plane_num' not in cfg:
        cfg['3D:from_plane_num'] = 0

    if '3D:to_plane_num' not in cfg:
        cfg['3D:to_plane_num'] = len(coo) - 1

    return spots, coo, cfg


if __name__ == "__main__":

    # set up the logger
    attach_to_log()

    # read 2D some demo data
    _iss_spots_2D = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'iss', 'spots.csv'))
    _coo_2D = load_npz(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'))

    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(float).astype(np.uint32)

    # read 3D some demo data
    _iss_spots_3D = pd.read_csv(r"E:\data\Anne\220308 50umCF seq atto425 DY520XL MS002\spots_yxz.csv")
    _iss_spots_3D = _iss_spots_3D.assign(z_stack=_iss_spots_3D.z)
    _iss_spots_3D = _iss_spots_3D[['y', 'x', 'z_stack', 'Gene']]

    _coo_3D = np.load(r"E:\data\Anne\220308 50umCF seq atto425 DY520XL MS002\masks_2D_stiched_fullsize.npz",  allow_pickle=True)['arr_0']
    _coo_3D = [coo_matrix(d) for d in _coo_3D]

    # main task
    opts_2D = {'save_data': True, 'nNeighbors': 3, 'MisreadDensity': 0.00001,'is_3D': False}
    opts_3D={'save_data': True,
             'Inefficiency': 0.0001,
             '3D:from_plane_num': 20,
             '3D:to_plane_num': 33,
             'MisreadDensity': 1e-05,
             'is_3D': True,
             'nNeighbors': 6,
          }

    fit(_iss_spots_2D, _coo_2D, scRNAseq=_scRNAseq, opts=opts_2D)
    # fit(_iss_spots_3D, _coo_3D, scRNAseq=_scRNAseq, opts=opts_3D)

