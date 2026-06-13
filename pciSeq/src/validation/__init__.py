"""
Simplified validation system for pciSeq.

This module provides a clean, single-function entry point for validating
all pciSeq inputs (spots, coo, scRNA data, configuration).

Usage:
    >>> from pciSeq.src.validation import validate_inputs
    >>> spots, coo, scdata, cfg = validate_inputs(spots, coo, scdata, opts)
"""
from typing import Tuple, Optional, Dict, Any
import pandas as pd
from scipy.sparse import coo_matrix
from .config import Config
from .validator import Validator


def validate_inputs(
    spots: pd.DataFrame,
    coo: Any,
    scdata: Optional[pd.DataFrame] = None,
    opts: Optional[Dict[str, Any]] = None
) -> Tuple[pd.DataFrame, list, Optional[pd.DataFrame], Dict]:
    """
    Validate all pciSeq inputs in one simple function call.

    This is the main entry point for the validation system. It:
        1. Creates and configures a Config object
        2. Sets runtime attributes (is3D)
        3. Runs comprehensive validation via Validator
        4. Returns validated, cleaned inputs

    Args:
        spots: DataFrame with spot data. Required columns: gene_name, x, y
               For 3D data, also needs z_plane column.
        coo: Segmentation label image(s). Either:
             - Single scipy.sparse.coo_matrix (for 2D)
             - List of coo_matrices (for 3D, one per z-plane)
        scdata: Optional single-cell RNA-seq reference data.
                Should be DataFrame with genes as index, cell types as columns.
        opts: Optional configuration overrides.
              Any key in pciSeq.config.DEFAULT can be overridden.

    Returns:
        tuple: (spots, coo, scdata, config_dict)
            - spots: Validated, cleaned DataFrame
            - coo: List of coo_matrices
            - scdata: Validated scRNA DataFrame or None
            - config_dict: Complete configuration dictionary

    Raises:
        TypeError: If inputs have incorrect types
        ValueError: If inputs have invalid values or structure
    """
    # Step 1: Create and configure Config
    config = Config(opts)
    config.set_runtime_attrs(coo)

    # Step 2: Run comprehensive validation
    validator = Validator(spots, coo, scdata, dict(config))
    return validator.validate_all()


# Convenience imports
__all__ = ['validate_inputs']