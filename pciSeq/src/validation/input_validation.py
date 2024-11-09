from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Union, List
from typing import get_origin, get_args
from .config_manager import ConfigManager
import pandas as pd
import numpy as np
from pandas import DataFrame
from scipy.sparse import coo_matrix
import logging

validation_logger = logging.getLogger(__name__)


@dataclass
class ValidatedInputs:
    """Container for validated inputs"""
    spots: pd.DataFrame
    coo: coo_matrix
    scdata: Optional[pd.DataFrame]
    config: ConfigManager


def validate_inputs(
        spots: pd.DataFrame,
        coo: coo_matrix,
        scdata: Optional[pd.DataFrame],
        config: ConfigManager
) -> tuple[DataFrame, coo_matrix, Union[DataFrame, None], Dict]:
    """Validate all inputs and return validated versions"""

    # Validate spots
    _validate_spots(spots)
    spots = _process_spots(spots.copy(), scdata, config)

    # Validate coo matrix
    _validate_coo(coo)

    # Validate single cell data if present
    if scdata is not None:
        _validate_scdata(scdata)

    # Validate config
    config = _process_config(config)

    out = ValidatedInputs(spots=spots, coo=coo, scdata=scdata, config=config)

    return out.spots, out.coo, out.scdata, out.config.to_dict()


def _validate_spots(spots: pd.DataFrame) -> None:
    """Validate spots dataframe structure"""
    if not isinstance(spots, pd.DataFrame):
        raise TypeError("Spots should be passed-in as a dataframe")

    if set(spots.columns) != {'Gene', 'x', 'y'}:
        raise ValueError("Spots dataframe must have columns ['Gene', 'x', 'y']")


def _validate_coo(coo: coo_matrix) -> None:
    """Validate coo matrix"""
    if not isinstance(coo, coo_matrix):
        raise TypeError('The segmentation masks should be passed-in as a coo_matrix')


def _validate_scdata(scdata: pd.DataFrame) -> None:
    """Validate single cell data"""
    if not isinstance(scdata, pd.DataFrame):
        raise TypeError("Single cell data should be passed-in as a dataframe")


def _process_spots(spots: pd.DataFrame, scdata: Optional[pd.DataFrame], config: 'ConfigManager') -> pd.DataFrame:
    """Process and validate spots data"""
    # Cast datatypes
    spots = spots.astype({
        'Gene': str,
        'x': np.float32,
        'y': np.float32
    })

    # Remove genes not in single cell data if present
    if scdata is not None:
        if not set(spots.Gene).issubset(scdata.index):
            spots = _purge_spots(spots, scdata)

    return spots


def _purge_spots(spots: pd.DataFrame, scdata: pd.DataFrame) -> pd.DataFrame:
    """Remove spots with genes not found in single cell data"""
    drop_spots = list(set(spots.Gene) - set(scdata.index))
    validation_logger.warning(f'Found {len(drop_spots)} genes that are not included in the single cell data')
    idx = ~np.in1d(spots.Gene, drop_spots)
    spots = spots.iloc[idx]
    validation_logger.warning(f'Removed from spots: {drop_spots}')
    return spots


def _process_config(config: 'ConfigManager') -> 'ConfigManager':
    """Process and validate configuration"""
    type_validations = {
        'exclude_genes': List[str],
        'max_iter': int,
        'CellCallTolerance': float,
        'rGene': int,
        'Inefficiency': float,
        'InsideCellBonus': Union[bool, int, float],
        'MisreadDensity': Union[float, Dict[str, float]],
        'SpotReg': float,
        'nNeighbors': int,
        'rSpot': int,
        'save_data': bool,
        'output_path': str,
        'launch_viewer': bool,
        'launch_diagnostics': bool,
        'is_redis_running': bool,
        'cell_radius': Optional[float],
        'cell_type_prior': str,
        'mean_gene_counts_per_class': int,
        'mean_gene_counts_per_cell': int
    }

    # check whether values have the expected type.
    # Note. This is not 100% safe. I will check the top-level type,
    # For 'exclude_genes' which ia a list of strings, it will
    # only check it is a list
    for attr_name, expected_type in type_validations.items():
        value = getattr(config, attr_name)
        origin_type = get_origin(expected_type)
        validation_logger.info(f'Doing {attr_name}')

        if origin_type is Union:
            allowed_types = get_args(expected_type)
            validation_logger.info(f'Allowed types {allowed_types}')
            allowed_types = tuple([get_origin(d) if get_origin(d) else d for d in allowed_types])
            if not isinstance(value, allowed_types):
                raise TypeError(f"'{attr_name}' must be one of {allowed_types}, got {type(value)}")
        else:
            origin_type = origin_type or expected_type
            validation_logger.info(f'Allowed types {origin_type}')
            if not isinstance(value, origin_type):
                raise TypeError(f"'{attr_name}' must be of type {expected_type}, got {type(value)}")

    if not isinstance(config.InsideCellBonus, (bool, int, float)):
        raise TypeError("'InsideCellBonus' must be a boolean, integer, or float")

    if not isinstance(config.MisreadDensity, (dict, int, float)):
        raise TypeError("'MisreadDensity' must be a dictionary, integer, or float")

    # Validate cell_type_prior
    if config.cell_type_prior.lower() not in ['uniform', 'weighted']:
        raise ValueError("'cell_type_prior' should be either 'uniform' or 'weighted'")

    # Convert to lowercase
    config.cell_type_prior = config.cell_type_prior.lower()

    # Process InsideCellBonus
    if config.InsideCellBonus is True:
        config.InsideCellBonus = 2
        validation_logger.warning('InsideCellBonus was passed-in as True. Overriding with the default value of 2')

    # Validate MisreadDensity
    if isinstance(config.MisreadDensity, dict):
        if 'default' not in config.MisreadDensity:
            raise ValueError("When MisreadDensity is a dictionary, it must contain a 'default' key")
    elif isinstance(config.MisreadDensity, (int, float)):
        config.MisreadDensity = {'default': config.MisreadDensity}
    else:
        raise ValueError("MisreadDensity must be either a number or a dictionary with a 'default' key")

    return config