"""
Unified input validation for pciSeq.

This module provides a single Validator class that handles all input
validation in a clear, sequential pipeline.
"""
from typing import Tuple, Optional, Dict, Any, List
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import logging

logger = logging.getLogger(__name__)


class Validator:
    """
    Single class for validating all pciSeq inputs.

    Handles validation of spots, coo matrices, scRNA data, and configuration
    in a clear, sequential pipeline with comprehensive error messages.

    Example:
        >>> validator = Validator(spots, coo, scdata, config)
        >>> spots, coo, scdata, cfg = validator.validate_all()
    """

    def __init__(
        self,
        spots: pd.DataFrame,
        coo: Any,
        scdata: Optional[pd.DataFrame],
        config: Dict[str, Any]
    ):
        """
        Initialize validator with all inputs.

        Args:
            spots: DataFrame with spot data
            coo: COO matrix or list of matrices (segmentation)
            scdata: Optional scRNA reference DataFrame
            config: Configuration dictionary
        """
        self.spots = spots
        self.coo = coo
        self.scdata = scdata
        self.config = config

    def validate_all(self) -> Tuple[pd.DataFrame, List[coo_matrix], Optional[pd.DataFrame], Dict]:
        """
        Run all validation steps in sequence.

        Validation pipeline:
            1. Type validation & input normalization
            2. Data structure validation
            3. Data cleaning & processing
            4. Configuration validation

        Returns:
            tuple: (validated_spots, validated_coo, validated_scdata, validated_config)

        Raises:
            TypeError: If inputs have wrong types
            ValueError: If inputs have invalid values or structure
        """
        logger.info("Starting validation pipeline")

        # Step 1: Validate types and normalize formats
        self._validate_types()
        self._normalize_inputs()

        # Step 2: Validate data structures
        self._validate_spots_schema()
        self._validate_coo_structure()
        if self.scdata is not None:
            self._validate_scdata_schema()

        # Step 3: Clean and process data
        self._clean_spots()
        self._reconcile_genes()

        # Step 4: Validate configuration
        self._validate_config()

        logger.info("Validation complete")
        return self.spots, self.coo, self.scdata, self.config

    # ============================================================================
    # STEP 1: Type validation & normalization
    # ============================================================================

    def _validate_types(self) -> None:
        """
        Check that inputs have correct basic types.

        Raises:
            TypeError: If spots or scdata are not DataFrames
        """
        if not isinstance(self.spots, pd.DataFrame):
            raise TypeError("spots must be a pandas DataFrame")

        if self.scdata is not None and not isinstance(self.scdata, pd.DataFrame):
            raise TypeError("scdata must be a pandas DataFrame")

    def _normalize_inputs(self) -> None:
        """
        Normalize inputs to expected internal formats.

        - Converts single coo_matrix to list format
        - Validates that all coo elements are matrices

        Raises:
            TypeError: If coo has invalid structure
        """
        # Normalize coo to list format
        if isinstance(self.coo, coo_matrix):
            self.coo = [self.coo]

        if not isinstance(self.coo, list):
            raise TypeError("coo must be a coo_matrix or list of coo_matrix")

        if not all(isinstance(m, coo_matrix) for m in self.coo):
            raise TypeError("All coo elements must be coo_matrix instances")

    # ============================================================================
    # STEP 2: Data structure validation
    # ============================================================================

    def _validate_spots_schema(self) -> None:
        """
        Validate spots DataFrame has required columns and handle 2D/3D differences.

        For 2D data:
            - Renames 'Gene' to 'gene_name' if needed
            - Adds z_plane=0 if missing

        Raises:
            ValueError: If required columns are missing
        """
        # Handle 2D data conventions
        if not self.config['is3D']:
            # Legacy column name support
            if 'Gene' in self.spots.columns:
                self.spots = self.spots.rename(columns={'Gene': 'gene_name'})

            # Add z_plane for 2D data
            if 'z_plane' not in self.spots.columns:
                self.spots['z_plane'] = 0

        # Check required columns
        required = {'gene_name', 'x', 'y', 'z_plane'}
        missing = required - set(self.spots.columns)

        if missing:
            raise ValueError(
                f"spots DataFrame missing required columns: {missing}. "
                f"Found columns: {list(self.spots.columns)}"
            )

    def _validate_coo_structure(self) -> None:
        """
        Validate coo matrices structure.

        The original validation only checked that coo was a non-empty list.
        Empty matrices within the list were allowed (common for edge planes in 3D).
        Logs warnings for empty planes.
        """
        if len(self.coo) == 0:
            raise ValueError("coo cannot be an empty list")

        # Check for empty planes and warn (but don't fail)
        empty_planes = []
        for i, matrix in enumerate(self.coo):
            if matrix.nnz == 0:  # nnz = number of non-zero elements
                empty_planes.append(i)

        if empty_planes:
            logger.warning(
                f"Found {len(empty_planes)} empty plane(s) with no cells: {empty_planes}"
            )

    def _validate_scdata_schema(self) -> None:
        """
        Validate single cell reference data is not empty.

        Raises:
            ValueError: If scdata is empty
        """
        if self.scdata.empty:
            raise ValueError("scdata cannot be empty")

    # ============================================================================
    # STEP 3: Data cleaning & processing
    # ============================================================================

    def _clean_spots(self) -> None:
        """
        Clean and typecast spots data to correct dtypes.

        Ensures:
            - gene_name: str
            - x, y, z_plane: float32
        """
        self.spots = self.spots.astype({
            'gene_name': str,
            'x': np.float32,
            'y': np.float32,
            'z_plane': np.float32
        })

    def _reconcile_genes(self) -> None:
        """
        Remove spots with genes not present in scRNA reference.

        Only runs if scRNA data is provided. Logs detailed information
        about removed genes and remaining spots.
        """
        if self.scdata is None:
            return

        spot_genes = set(self.spots['gene_name'])
        scdata_genes = set(self.scdata.index)

        # Check if all spot genes are in scRNA data
        if spot_genes.issubset(scdata_genes):
            logger.info("All spot genes found in scRNA reference")
            return

        # Find missing genes
        missing_genes = spot_genes - scdata_genes
        logger.warning(
            f"Found {len(missing_genes)} genes in spots but not in scRNA reference"
        )
        logger.warning(f"Missing genes: {sorted(missing_genes)}")

        # Filter spots
        initial_count = len(self.spots)
        mask = self.spots['gene_name'].isin(scdata_genes)
        self.spots = self.spots[mask].reset_index(drop=True)

        removed_count = initial_count - len(self.spots)
        logger.info(
            f"Removed {removed_count} spots ({removed_count/initial_count*100:.1f}%). "
            f"Kept {len(self.spots)} spots"
        )

    # ============================================================================
    # STEP 4: Config validation
    # ============================================================================

    def _validate_config(self) -> None:
        """
        Validate and normalize configuration values.

        Handles:
            - Rejection of deprecated parameters (exclude_planes)
            - cell_type_prior validation and normalization
            - InsideCellBonus conversion (True to 2)
            - Dict parameter normalization (MisreadDensity, priors)
        """
        cfg = self.config

        # Check for deprecated/unsupported parameters
        if 'exclude_planes' in cfg and cfg['exclude_planes'] is not None:
            raise ValueError(
                "The 'exclude_planes' parameter is no longer supported due to coordinate system "
                "complexity. Please filter your input data (spots and coo matrices) before "
                "passing them to pciSeq. Use spatial filtering (x, y, z) on your input data instead."
            )

        # Validate and normalize cell_type_prior
        if cfg['cell_type_prior'].lower() not in ['uniform', 'weighted']:
            raise ValueError(
                "cell_type_prior must be 'uniform' or 'weighted', "
                f"got '{cfg['cell_type_prior']}'"
            )
        cfg['cell_type_prior'] = cfg['cell_type_prior'].lower()

        # Normalize InsideCellBonus boolean to numeric
        if cfg['InsideCellBonus'] is True:
            cfg['InsideCellBonus'] = 2
            logger.warning("InsideCellBonus=True converted to default value of 2")

        # Normalize dict parameters
        cfg['MisreadDensity'] = self._ensure_dict(
            cfg['MisreadDensity'], 'MisreadDensity'
        )
        cfg['cell_centroid_prior'] = self._ensure_dict(
            cfg['cell_centroid_prior'], 'cell_centroid_prior'
        )
        cfg['cell_cov_prior'] = self._ensure_dict(
            cfg['cell_cov_prior'], 'cell_cov_prior'
        )

    @staticmethod
    def _ensure_dict(param: Any, name: str) -> Dict[str, Any]:
        """
        Convert scalar to dict or validate dict has 'default' key.

        This allows users to pass either:
            - A scalar: 0.1 → {'default': 0.1}
            - A dict: {'default': 0.1, 'Plp1': 0.2}

        Args:
            param: Parameter value (scalar or dict)
            name: Parameter name (for error messages)

        Returns:
            dict: Dictionary with 'default' key

        Raises:
            ValueError: If param is invalid type or dict missing 'default'
        """
        # Convert scalar to dict
        if isinstance(param, (int, float)):
            return {'default': param}

        # Validate dict has 'default' key
        if isinstance(param, dict):
            if 'default' not in param:
                raise ValueError(
                    f"{name} dictionary must contain 'default' key. "
                    f"Found keys: {list(param.keys())}"
                )
            return param

        # Invalid type
        raise ValueError(
            f"{name} must be a number or dict with 'default' key, "
            f"got {type(param).__name__}"
        )