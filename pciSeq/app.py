import pandas as pd
from typing import Tuple, Optional, Dict, Any
from .src.validation.config_manager import ConfigManager
from .src.validation.input_validation import validate_inputs
from .src.core.main import VarBayes
from .src.core.utils import write_data, pre_launch
from .src.viewer.run_flask import flask_app_start
from .src.preprocess.spot_labels import stage_data
from .src.diagnostics.launch_diagnostics import launch_dashboard
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
    try:
        # 1. parse/check the arguments
        spots, coo, scRNAseq, opts = parse_args(*args, **kwargs)

        # 2. Create and validate config
        cfg_man = ConfigManager.from_opts(opts)
        spots, coo, scdata, cfg = validate_inputs(spots, coo, scRNAseq, cfg_man)

        # 3. Use validated inputs and prepare the data
        app_logger.info('Preprocessing data')
        _cells, cellBoundaries, _spots = stage_data(spots, coo)

        # 4. launch the diagnostics
        if cfg['launch_diagnostics'] and cfg.get('is_redis_running', False):  # Fixed typo in condition
            app_logger.info('Launching the diagnostics dashboard')
            launch_dashboard()

        # 5. cell typing
        cellData, geneData, varBayes = cell_type(_cells, _spots, scdata, cfg)

        # 6. Save data and launch viewer if needed
        if (cfg['save_data'] and varBayes.has_converged) or cfg['launch_viewer']:
            write_data(cellData, geneData, cellBoundaries, varBayes, cfg)

            if cfg['launch_viewer']:
                dst = pre_launch(cellData, coo, scRNAseq, cfg)
                flask_app_start(dst)

        app_logger.info('Done')
        return cellData, geneData

    except Exception as e:
        app_logger.error(f"Error in fit function: {str(e)}")
        raise


def cell_type(
        cells: pd.DataFrame,
        spots: pd.DataFrame,
        scRNAseq: Optional[pd.DataFrame],
        config: Dict[str, Any]
) -> Tuple[pd.DataFrame, pd.DataFrame, VarBayes]:
    """
    Perform cell typing using Variational Bayes algorithm.

    Parameters
    ----------
    cells : pd.DataFrame
        Preprocessed cell data containing cell locations and boundaries
    spots : pd.DataFrame
        Preprocessed spot data containing gene expressions and coordinates
    scRNAseq : Optional[pd.DataFrame]
        Single-cell RNA sequencing reference data. Can be None if not using reference data
    config : Dict[str, Any]
        Configuration dictionary containing algorithm parameters

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, VarBayes]
        - cellData: DataFrame containing cell typing results
        - geneData: DataFrame containing gene assignment results
        - varBayes: The fitted VarBayes model instance

    Raises
    ------
    ValueError
        If input data is invalid or incompatible
    RuntimeError
        If cell typing algorithm fails to converge
    """
    try:
        app_logger.info('Initializing VarBayes model')
        varBayes = VarBayes(cells, spots, scRNAseq, config)

        app_logger.info('Starting cell typing algorithm')
        cellData, geneData = varBayes.run()

        if not varBayes.has_converged:
            app_logger.warning('Cell typing algorithm did not fully converge')

        return cellData, geneData, varBayes

    except Exception as e:
        app_logger.error(f"Error during cell typing: {str(e)}")
        raise RuntimeError(f"Cell typing failed: {str(e)}") from e


def parse_args(*args, **kwargs) -> Tuple[pd.DataFrame, Any, Optional[pd.DataFrame], Optional[Dict]]:
    """Parse and validate input arguments.

    Returns:
        Tuple containing (spots, coo, scRNAseq, opts)

    Raises:
        ValueError: If required arguments are missing or invalid
    """
    # Check if we have the minimum required arguments
    if not {'spots', 'coo'}.issubset(set(kwargs)) and len(args) < 2:
        raise ValueError('Need to provide the spots and the coo matrix either as keyword '
                         'arguments or as the first and second positional arguments.')

    # Get spots from kwargs if present, otherwise from args
    spots = kwargs['spots'] if 'spots' in kwargs else args[0]

    # Get coo from kwargs if present, otherwise from args
    coo = kwargs['coo'] if 'coo' in kwargs else args[1]

    # Optional arguments
    scRNAseq = kwargs.get('scRNAseq', None)
    opts = kwargs.get('opts', None)

    if spots is None or coo is None:
        raise ValueError("Both 'spots' and 'coo' must be provided")

    return spots, coo, scRNAseq, opts






