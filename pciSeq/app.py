import pandas as pd
from scipy.sparse import coo_matrix
import numpy as np
from typing import Tuple, Optional, Dict, Any
from .src.validation.config_manager import ConfigManager
from .src.validation.input_validation import InputValidator
from .src.core.main import VarBayes
from .src.core.utils.cell_utils import recover_original_labels
from .src.core.utils.io_utils import write_data
from .src.viewer.utils import pre_launch
from .src.viewer.run_flask import flask_app_start
from .src.preprocess.main import stage_data
import logging

app_logger = logging.getLogger(__name__)


def fit(*args, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for pciSeq cell typing analysis.

    Parameters
    ----------
    *args : tuple
        Positional arguments:
        - args[0]: pd.DataFrame containing spot data
        - args[1]: scipy.sparse.coo_matrix containing gene expression data

    **kwargs : dict
        Keyword arguments (preferred method):
        spots : pd.DataFrame
            Spot data with required columns:
            - 'gene_name': Name of the gene
            - 'x': X coordinate
            - 'y': Y coordinate
            - 'z_plane': Z plane (optional for 3D data)

        coo : List[scipy.sparse.coo_matrix]
            List of sparse matrices containing gene expression data.
            Length > 1 indicates 3D data

        scRNAseq : pd.DataFrame, optional
            Single-cell RNA sequencing reference data.
            Used for cell type annotation if provided

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - cellData: DataFrame containing cell typing results and metadata
        - geneData: DataFrame containing gene assignment results

    Raises
    ------
    ValueError
        If required arguments (spots and coo) are missing or invalid
    RuntimeError
        If cell typing algorithm fails to converge

    Notes
    -----
    The function can be called either with positional arguments (spots, coo)
    or with keyword arguments. If both are provided, keyword arguments take precedence.
    """
    try:
        # 1. parse/check the arguments
        spots, coo, scRNAseq, opts = parse_args(*args, **kwargs)

        # 2. Create and validate config
        cfg_man = ConfigManager.from_opts(opts)
        cfg_man.set_runtime_attributes(coo)
        spots, coo, scdata, cfg = InputValidator.validate(spots, coo, scRNAseq, cfg_man)

        # 3. Use validated inputs and prepare the data
        app_logger.info('Preprocessing data')
        _cells, cellBoundaries, _spots, label_map = stage_data(spots, coo, cfg)
        cfg['remapping'] = label_map

        # 5. cell typing (diagnostics are now handled inside VarBayes)
        cellData, geneData, varBayes = cell_type(_cells, _spots, scdata, cfg)

        # 6 if labels have been remapped, switch to the original ones
        if label_map is not None:
            cellData, geneData, cellBoundaries = recover_original_labels(cellData, geneData, cellBoundaries, label_map)

        # 7. Save data and launch viewer if needed
        if cfg['save_data'] or cfg['launch_viewer']:
            write_data(cellData, geneData, cellBoundaries, varBayes, cfg)

            if cfg['launch_viewer']:
                dst = pre_launch(cellData, geneData, coo, scRNAseq, cfg)
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
    if isinstance(coo, np.ndarray):
        coo = [coo_matrix(d) for d in coo]

    # Optional arguments
    scRNAseq = kwargs.get('scRNAseq', None)
    opts = kwargs.get('opts', None)

    if spots is None or coo is None:
        raise ValueError("Both 'spots' and 'coo' must be provided")

    return spots, coo, scRNAseq, opts






