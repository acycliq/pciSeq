import pandas as pd
from typing import Tuple
from .src.core.config_manager import ConfigManager
from .src.core.validation import validate_inputs
from pciSeq.src.core.main import VarBayes
from pciSeq.src.core.utils import write_data, pre_launch
from pciSeq.src.viewer.run_flask import flask_app_start
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.diagnostics.launch_diagnostics import launch_dashboard
import logging

app_logger = logging.getLogger(__name__)


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
            Name: Genenames, dtype: Object, array-like of the genes assigned to the cell
            Name: CellGeneCount, dtype: Object,array-like of the corresponding gene counts
            Name: ClassName, dtype: Object, array-like of the genes probable classes for the cell
            Name: Prob, dtype: Object, array-like of the corresponding cell class probabilities

    geneData : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Gene, dtype: string, The gene name.
            Name: Gene_id, dtype: int64, The gene id, the position of the gene if all genes are sorted.
            Name: x, dtype: int64, X-axis coordinate of the spot
            Name: y, dtype: int64, Y-axis coordinate of the spot
            Name: neighbour, dtype: int64, the label of the cell which is more likely to 'raise' the spot. If zero then the spot is a misread.
            Name: neighbour_array, dtype: Object, array-like with the labels of the 4 nearest cell. The last is always the background and has label=0
            Name: neighbour_prob, dtype: Object, array-like with the prob the corresponding cell from neighbour_array has risen the spot.
    """

    # 1. parse/check the arguments
    spots, coo, scRNAseq, opts = parse_args(*args, **kwargs)

    # 2. Create config
    cfg_man = ConfigManager.from_opts(opts)

    # 3. validate inputs
    spots, coo, scdata, cfg = validate_inputs(spots, coo, scRNAseq, cfg_man)

    # 4. Use validated inputs and prepare the data
    app_logger.info('Preprocessing data')
    _cells, cellBoundaries, _spots = stage_data(spots, coo)

    # 5. launch the diagnostics
    if cfg['launch_diagnostics'] and ['cfg.is_redis_running']:
        app_logger.info('Launching the diagnostics dashboard')
        launch_dashboard()

    # 6. cell typing
    cellData, geneData, varBayes = cell_type(_cells, _spots, scdata, cfg)

    # 7. save to the filesystem
    if (cfg['save_data'] and varBayes.has_converged) or cfg['launch_viewer']:
        write_data(cellData, geneData, cellBoundaries, varBayes, cfg)

    # 8. do the viewer if needed
    if cfg['launch_viewer']:
        dst = pre_launch(cellData, coo, scRNAseq, cfg)
        flask_app_start(dst)

    app_logger.info('Done')
    return cellData, geneData


def cell_type(_cells, _spots, scRNAseq, ini):
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)

    app_logger.info('Start cell typing')
    cellData, geneData = varBayes.run()
    return cellData, geneData, varBayes


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
                                   'and second args to the fit() method.'
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






