from dataclasses import dataclass
from typing import Dict, Optional, Union
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

    # Validate and process config
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