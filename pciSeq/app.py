import pandas as pd
from scipy.sparse import coo_matrix
import numpy as np
from typing import Tuple, Optional, Dict, Any
from .src.validation import validate_inputs
from .src.core.main import VarBayes
from .src.core.utils.cell_utils import recover_original_labels
from .src.core.io import write_data
from .src.preprocess.main import stage_data
import logging

logger = logging.getLogger(__name__)


def fit(*args, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for pciSeq cell typing analysis.

    ``spots`` and ``coo`` may be given as the first two positional arguments or as
    keywords. Keyword form is preferred.

    Parameters
    ----------
    spots : pd.DataFrame
        The spots to assign. Needs the columns 'gene_name', 'x' and 'y', plus
        'z_plane' for 3D data. Optional 'score' and 'intensity' columns from the
        spot caller are carried through to geneData; they default to 1.0 and
        nothing in the model reads them.

    coo : list of scipy.sparse.coo_matrix, or numpy array
        The label image: each pixel holds the label of the cell it belongs to, and
        background pixels are labelled 0. Pass it as one sparse matrix per z-plane,
        or as a dense numpy array of shape (height, width) or (planes, height,
        width), which is converted to sparse internally. A single coo_matrix is also
        accepted. More than one plane is treated as 3D. Sparse matrices are modified
        in place: labels that are not sequential are renumbered, and with
        `remove_flat_cells` the cells on a single plane are zeroed; pass a copy to
        keep the original. A dense array is left untouched.

    scRNAseq : pd.DataFrame
        Cell type definitions: mean expression per gene and cell type, genes as
        rows and cell types as columns. Required.

    opts : dict, optional
        Config values to override, e.g. {'max_iter': 500}. See the
        configuration page for the full list of keys and their defaults.

    Returns
    -------
    cellData : pd.DataFrame
        One row per cell: the class probabilities, the gene counts and the cell
        geometry.

    geneData : pd.DataFrame
        One row per spot: which cell it was assigned to and with what probability.

    Raises
    ------
    ValueError
        If spots, coo or scRNAseq are missing or malformed, or an option value is
        out of range.
    KeyError
        If opts has a key that is not a config option.
    TypeError
        If an input or an option value has the wrong type.
    Exception
        Any error raised inside the model is logged and re-raised unchanged. Not
        converging is not an error: the loop runs to max_iter, logs the convergence
        status and returns its results.

    Notes
    -----
    `spots` and `coo` are given either both as keywords or as the first two positional
    arguments. `scRNAseq` and `opts` are keyword only.
    """
    viewer = None  # Track realtime viewer for cleanup
    try:
        # 1. parse/check the arguments
        spots, coo, scRNAseq, opts = parse_args(*args, **kwargs)

        # 2. Validate all inputs (spots, coo, scRNA, config)
        spots, coo, scdata, cfg = validate_inputs(spots, coo, scRNAseq, opts)

        # 3. Start realtime viewer if requested
        viewer = realtime_viewer_ini(cfg)

        # 4. Use validated inputs and prepare the data
        logger.info('Preprocessing data')
        _cells, borders_future, _spots, label_map = stage_data(spots, coo, cfg)

        # 5. cell typing (diagnostics are now handled inside VarBayes)
        cellData, geneData, varBayes = cell_type(_cells, _spots, scdata, cfg)

        # 6. Resolve borders (blocks only if extraction hasn't finished yet)
        cellBoundaries, cellBoundaries_list = borders_future.result()

        # 7 if labels have been remapped, switch to the original ones. It is a no-op
        # when there was no renumbering, ie when label_map is None.
        cellData, geneData, cellBoundaries, cellBoundaries_list = recover_original_labels(
            cellData, geneData, cellBoundaries, cellBoundaries_list, label_map)

        # 8. Save data
        if cfg['save_data']:
            write_data(cellData, geneData, cellBoundaries, cellBoundaries_list, coo, varBayes, cfg)

        logger.info('Done')
        return cellData, geneData

    except Exception as e:
        logger.error(f"Error in fit function: {str(e)}")
        raise
    finally:
        # Cleanup realtime viewer if it was started
        if viewer is not None:
            try:
                viewer.stop()
                logger.info('Stopped realtime viewer')
            except Exception as e:
                logger.warning(f'Failed to stop realtime viewer: {e}')


def cell_type(
        cells: pd.DataFrame,
        spots: pd.DataFrame,
        scRNAseq: pd.DataFrame,
        config: Dict[str, Any]
) -> Tuple[pd.DataFrame, pd.DataFrame, VarBayes]:
    """
    Perform cell typing using Variational Bayes algorithm.

    Parameters
    ----------
    cells : pd.DataFrame
        One row per cell, with columns 'label', 'area' and the centroid 'x0', 'y0',
        'z0'. stage_data produces it.
    spots : pd.DataFrame
        One row per spot, with columns 'x', 'y', 'z', 'plane_id', 'gene_name', 'score',
        'intensity' and 'label', the label of the cell the spot lies in, 0 for none.
        stage_data produces it.
    scRNAseq : pd.DataFrame
        Cell type definitions: mean expression per gene and cell type, genes as
        rows and cell types as columns. Required.
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
    Exception
        Any other error raised inside the model is logged and re-raised unchanged.
        Not converging is not an error, it only logs a warning.
    """
    try:
        # Extract callback from config BEFORE creating VarBayes
        # This prevents the callback from being stored in VarBayes.config
        callback = config.pop('realtime_viewer_callback', None)

        logger.info('Initializing VarBayes model')
        varBayes = VarBayes(cells, spots, scRNAseq, config)

        # Wire real-time viewer callback if provided
        if callback is not None:
            # If callback is a bound method (e.g., viewer.send_update),
            # set the VarBayes reference so it can access cells data
            if hasattr(callback, '__self__'):
                viewer_instance = callback.__self__
                viewer_instance._varbayes_ref = varBayes

            varBayes.on_iteration_callback = callback
            logger.info('Real-time viewer callback enabled')

        logger.info('Starting cell typing algorithm')
        cellData, geneData = varBayes.run()

        if not varBayes.has_converged:
            logger.warning('Cell typing algorithm did not fully converge')

        return cellData, geneData, varBayes

    except Exception as e:
        logger.error(f"Error during cell typing: {str(e)}")
        # re-raise it as it is. Wrapping everything in a RuntimeError threw the type
        # away, so bad input and a crash inside the model looked the same to the caller.
        raise


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
        # here coo is the dense label image. A 2D one is a single plane. Looping over
        # it would hand back one "plane" per pixel row, which is what used to happen
        # and nothing complained about it.
        if coo.ndim == 2:
            coo = [coo_matrix(coo)]
        elif coo.ndim == 3:
            coo = [coo_matrix(d) for d in coo]
        else:
            raise ValueError(f"coo as a numpy array has to be 2D (height, width) or "
                             f"3D (planes, height, width), got {coo.ndim}D")

    # Optional arguments
    scRNAseq = kwargs.get('scRNAseq', None)
    opts = kwargs.get('opts', None)

    if spots is None or coo is None:
        raise ValueError("Both 'spots' and 'coo' must be provided")

    return spots, coo, scRNAseq, opts


def realtime_viewer_ini(cfg):
    if cfg.get("realtime_viewer", False):
        from .src.realtime_viewer import RealtimeViewerServer

        port = cfg.get("realtime_viewer_port", 5001)
        fixed_radius = cfg.get("realtime_viewer_fixed_radius", None)

        viewer = RealtimeViewerServer(port=port, fixed_radius=fixed_radius)
        viewer.start()
        cfg["realtime_viewer_callback"] = viewer.send_update
        logger.info(f"Started realtime viewer on port {port}")
        # hand it back, fit() keeps it so the finally block can stop it. Without
        # this the viewer was built here and dropped, so fit's `viewer` stayed
        # None and the cleanup never ran.
        return viewer
    return None





