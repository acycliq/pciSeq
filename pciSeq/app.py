import pandas as pd
from scipy.sparse import coo_matrix
import numpy as np
from typing import Tuple, Optional, Dict, Any
from .src.validation import validate_inputs
from .src.core.main import VarBayes
from .src.core.utils.cell_utils import recover_original_labels
from .src.core.utils.io_utils import write_data
from .src.viewer.utils import pre_launch
from .src.viewer.run_flask import flask_app_start
from .src.preprocess.main import stage_data
import logging

logger = logging.getLogger(__name__)


def fit(*args, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for pciSeq cell typing analysis.

    Parameters
    ----------
    *args : tuple
        Positional arguments:
        - args[0]: pd.DataFrame containing spot data
        - args[1]: list of scipy.sparse.coo_matrix (one per z-plane), the label image

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
        If the cell typing algorithm fails (non-convergence only logs a warning)

    Notes
    -----
    The function can be called either with positional arguments (spots, coo)
    or with keyword arguments. If both are provided, keyword arguments take precedence.
    """
    viewer = None  # Track realtime viewer for cleanup
    try:
        # 1. parse/check the arguments
        spots, coo, scRNAseq, opts = parse_args(*args, **kwargs)

        # 2. Validate all inputs (spots, coo, scRNA, config)
        spots, coo, scdata, cfg = validate_inputs(spots, coo, scRNAseq, opts)

        # 3. Start realtime viewer if requested (keep the handle so we can stop it)
        viewer = realtime_viewer_ini(cfg)

        # 4. Use validated inputs and prepare the data
        logger.info('Preprocessing data')
        _cells, borders_future, _spots = stage_data(spots, coo, cfg)

        # 5. cell typing (diagnostics are now handled inside VarBayes)
        cellData, geneData, varBayes = cell_type(_cells, _spots, scdata, cfg, viewer)

        # 6. Resolve borders (blocks only if extraction hasn't finished yet)
        cellBoundaries, cellBoundaries_list = borders_future.result()

        # 7 if labels have been remapped, switch to the original ones
        # (stage_data put label_map in cfg; that's the single home now)
        label_map = cfg.get('label_map')
        if label_map is not None:
            cellData, geneData, cellBoundaries, cellBoundaries_list = recover_original_labels(cellData, geneData, cellBoundaries, cellBoundaries_list, label_map)

        # 8. Save data and launch viewer if needed
        if cfg['save_data'] or cfg['launch_viewer']:
            write_data(cellData, geneData, cellBoundaries, cellBoundaries_list, varBayes, cfg)

            if cfg['launch_viewer']:
                dst = pre_launch(cellData, geneData, coo, scRNAseq, cfg)
                flask_app_start(dst)

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
        scRNAseq: Optional[pd.DataFrame],
        config: Dict[str, Any],
        viewer: Optional[Any] = None
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
    viewer : optional
        A running RealtimeViewerServer to stream iterations to, or None.
        When given, it is bound to the model via viewer.attach(varBayes).

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
        If the cell typing algorithm fails (non-convergence only logs a warning)
    """
    try:
        logger.info('Initializing VarBayes model')
        varBayes = VarBayes(cells, spots, scRNAseq, config)

        # Wire the realtime viewer to the model, if one is running.
        if viewer is not None:
            viewer.attach(varBayes)
            logger.info('Real-time viewer callback enabled')

        logger.info('Starting cell typing algorithm')
        cellData, geneData = varBayes.run()

        if not varBayes.has_converged:
            logger.warning('Cell typing algorithm did not fully converge')

        return cellData, geneData, varBayes

    except Exception as e:
        logger.error(f"Error during cell typing: {str(e)}")
        raise RuntimeError(f"Cell typing failed: {str(e)}") from e


def parse_args(*args, **kwargs) -> Tuple[pd.DataFrame, Any, Optional[pd.DataFrame], Optional[Dict]]:
    """Parse and validate input arguments.

    Returns:
        Tuple containing (spots, coo, scRNAseq, opts)

    Raises:
        ValueError: If required arguments are missing or invalid
    """
    # spots and coo can come either as the first two positional args or as
    # keywords. Resolve each on its own (keyword wins) so that a missing one
    # gives a clear error below instead of an IndexError or a wrong slot.
    spots = kwargs.get('spots', args[0] if len(args) > 0 else None)
    coo = kwargs.get('coo', args[1] if len(args) > 1 else None)
    scRNAseq = kwargs.get('scRNAseq', None)
    opts = kwargs.get('opts', None)

    if spots is None or coo is None:
        raise ValueError('Need to provide both spots and coo, either as the first two '
                         'positional arguments or as the spots= and coo= keyword arguments.')

    # a 3D stack can arrive as a single ndarray; turn it into the list of sparse
    # matrices (one per z-plane) that the rest of the pipeline expects.
    if isinstance(coo, np.ndarray):
        coo = [coo_matrix(d) for d in coo]

    return spots, coo, scRNAseq, opts


def realtime_viewer_ini(cfg):
    """Start the realtime viewer if the config asks for it.

    Returns the running viewer so the caller can stop it later, or None when
    the viewer is switched off. cell_type wires it into the model via
    viewer.attach(varBayes).
    """
    if not cfg.get("realtime_viewer", False):
        return None

    from .src.realtime_viewer import RealtimeViewerServer

    port = cfg.get("realtime_viewer_port", 5001)
    max_cells = cfg.get("realtime_viewer_max_cells", None)
    fixed_radius = cfg.get("realtime_viewer_fixed_radius", None)

    viewer = RealtimeViewerServer(
        port=port, max_cells=max_cells, fixed_radius=fixed_radius
    )
    viewer.start()
    logger.info(f"Started realtime viewer on port {port}")
    return viewer





