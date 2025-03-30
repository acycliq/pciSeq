from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Union, List
from typing import get_origin, get_args, get_type_hints
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import logging
from .config_manager import ConfigManager

validation_logger = logging.getLogger(__name__)


@dataclass
class InputValidator:
    """Container for pciSeq analysis inputs with validation methods"""
    spots: pd.DataFrame
    coo: coo_matrix
    scdata: Optional[pd.DataFrame]
    config: ConfigManager

    @classmethod
    def validate(
            cls,
            spots: pd.DataFrame,
            coo: coo_matrix,
            scdata: Optional[pd.DataFrame],
            config: ConfigManager
    ) -> tuple[pd.DataFrame, coo_matrix, Union[pd.DataFrame, None], Dict]:
        """
        Validate all inputs for pciSeq analysis.

        Parameters
        ----------
        spots : pd.DataFrame
            Dataframe containing spot data with required columns
        coo : coo_matrix
            List of sparse matrices containing image segmentation labeled cells
        scdata : Optional[pd.DataFrame]
            Single-cell reference data, if available
        config : ConfigManager
            Configuration object with algorithm parameters

        Returns
        -------
        tuple
            Validated versions of (spots, coo, scdata, config_dict)

        Raises
        ------
        TypeError
            If inputs are of incorrect type
        ValueError
            If input data is invalid or incompatible
        """
        # Validate spots
        spots = cls._validate_spots(spots, config)
        spots = cls._process_spots(spots.copy(), scdata, config)

        # Validate coo matrix
        coo = cls._validate_coo(coo)

        # Validate single cell data if present
        if scdata is not None:
            cls._validate_scdata(scdata)

        # Validate config
        config = cls._validate_config(config)

        out = cls(spots=spots, coo=coo, scdata=scdata, config=config)
        return out.spots, out.coo, out.scdata, out.config.to_dict()

    @staticmethod
    def _validate_spots(spots: pd.DataFrame, config: Dict) -> pd.DataFrame:
        """
        Validate spots dataframe structure and content.

        Parameters
        ----------
        spots : pd.DataFrame
            Dataframe to validate

        Raises
        ------
        TypeError
            If spots is not a pandas DataFrame
        ValueError
            If required columns are missing
        """
        if not isinstance(spots, pd.DataFrame):
            raise TypeError("Spots should be passed-in as a dataframe")

        if not config.is3D:
            spots = spots.rename(columns={"Gene": "gene_name"})
            spots = spots.assign(z_plane=0)

        required_columns = {'gene_name', 'x', 'y', 'z_plane'}
        if not required_columns.issubset(spots.columns):
            raise ValueError(f"Spots dataframe must have columns {required_columns}")
        return spots

    @staticmethod
    def _validate_coo(coo: coo_matrix) -> None:
        """
        Validate coo matrix format and content.

        Parameters
        ----------
        coo : List[coo_matrix]
            List of sparse matrices to validate

        Raises
        ------
        TypeError
            If coo is not a list of sparse matrices
        ValueError
            If matrices are empty or invalid
        """
        if isinstance(coo, coo_matrix):
            coo = [coo]
        if not isinstance(coo, list) or not all(isinstance(item, coo_matrix) for item in coo):
            raise TypeError('The segmentation masks should be passed-in as a List[coo_matrix]')
        return coo

    @staticmethod
    def _validate_scdata(scdata: pd.DataFrame) -> None:
        """
        Validate single cell reference data.

        Parameters
        ----------
        scdata : pd.DataFrame
            Single cell reference data to validate

        Raises
        ------
        TypeError
            If scdata is not a pandas DataFrame
        ValueError
            If data is empty or invalid
        """
        if not isinstance(scdata, pd.DataFrame):
            raise TypeError("Single cell data should be passed-in as a dataframe")

    @staticmethod
    def _process_spots(spots: pd.DataFrame, scdata: Optional[pd.DataFrame], config: ConfigManager) -> pd.DataFrame:
        """
        Process and clean spots data.

        Parameters
        ----------
        spots : pd.DataFrame
            Spots data to process
        scdata : Optional[pd.DataFrame]
            Single cell reference data if available
        config : ConfigManager
            Configuration object

        Returns
        -------
        pd.DataFrame
            Processed spots data
        """
        # Cast datatypes
        spots = spots.astype({
            'gene_name': str,
            'x': np.float32,
            'y': np.float32,
            'z_plane': np.float32
        })

        # Remove genes not in single cell data if present
        if scdata is not None:
            if not set(spots.gene_name).issubset(scdata.index):
                spots = InputValidator._purge_spots(spots, scdata)

        return spots

    @staticmethod
    def _purge_spots(spots: pd.DataFrame, scdata: pd.DataFrame) -> pd.DataFrame:
        """
        Remove spots with genes not found in single cell data.

        Parameters
        ----------
        spots : pd.DataFrame
            Spots data to purge
        scdata : pd.DataFrame
            Reference single cell data

        Returns
        -------
        pd.DataFrame
            Purged spots data
        """
        drop_spots = list(set(spots.gene_name) - set(scdata.index))
        validation_logger.warning(f'Found {len(drop_spots)} genes that are not included in the single cell data')
        idx = ~np.in1d(spots.gene_name, drop_spots)
        spots = spots.iloc[idx]
        validation_logger.warning(f'Removed from spots: {drop_spots}')
        return spots

    @staticmethod
    def _validate_config(config: ConfigManager) -> ConfigManager:
        """
        Validate configuration settings.

        Parameters
        ----------
        config : ConfigManager
            Configuration to validate

        Returns
        -------
        ConfigManager
            Validated configuration

        Raises
        ------
        TypeError
            If config values have incorrect types
        ValueError
            If config values are invalid
        """
        # Type validation
        type_validations = get_type_hints(ConfigManager)
        for attr_name, expected_type in type_validations.items():
            value = getattr(config, attr_name)
            origin_type = get_origin(expected_type)

            if origin_type is Union:
                allowed_types = get_args(expected_type)
                allowed_types = tuple([get_origin(d) if get_origin(d) else d for d in allowed_types])
                if not isinstance(value, allowed_types):
                    raise TypeError(f"'{attr_name}' must be one of {allowed_types}, got {type(value)}")
            else:
                origin_type = origin_type or expected_type
                if not isinstance(value, origin_type):
                    raise TypeError(f"'{attr_name}' must be of type {expected_type}, got {type(value)}")

        # Value validation
        # if config.is3D:
        #     config.InsideCellBonus = False
        #     validation_logger.warning('InsideCellBonus set to False for 3D data')

        if config.cell_type_prior.lower() not in ['uniform', 'weighted']:
            raise ValueError("'cell_type_prior' should be either 'uniform' or 'weighted'")

        config.cell_type_prior = config.cell_type_prior.lower()

        if config.InsideCellBonus is True:
            config.InsideCellBonus = 2
            validation_logger.warning('InsideCellBonus was passed-in as True. Overriding with default value of 2')

        config.MisreadDensity = InputValidator._dict_checker(config.MisreadDensity, "MisreadDensity")
        config.cell_centroid_prior = InputValidator._dict_checker(config.cell_centroid_prior, "cell_centroid_prior")
        config.cell_cov_prior = InputValidator._dict_checker(config.cell_cov_prior, "cell_cov_prior")

        return config

    @staticmethod
    def _dict_checker(param, param_name):
        """ Simple helper method """
        if isinstance(param, dict):
            if 'default' not in param:
                raise ValueError(f"{param_name} dictionary must contain a 'default' key")
        elif isinstance(param, (int, float)):
            param = {'default': param}
        else:
            raise ValueError(f"{param_name} must be either a number or a dictionary with a 'default' key")
        return param

