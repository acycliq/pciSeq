"""
configuration management for pciSeq.
"""

from typing import Dict, Any, Optional
from pciSeq import config
import logging

logger = logging.getLogger(__name__)


# Type specifications for validation
_TYPE_SPECS = {
    "exclude_genes": list,
    "max_iter": int,
    "CellCallTolerance": float,
    "rGene": int,
    "Inefficiency": float,
    "InsideCellBonus": (bool, int, float),
    "MisreadDensity": (float, dict),
    "cell_centroid_prior": (int, float, dict),
    "cell_cov_prior": (int, float, dict),
    "SpotReg": float,
    "nNeighbors": int,
    "rSpot": (int, float),
    "save_data": bool,
    "output_path": str,
    "launch_viewer": (bool, str),
    "launch_diagnostics": bool,
    "is_redis_running": bool,
    "cell_radius": (type(None), float),
    "cell_type_prior": str,
    "voxel_size": list,
    "exclude_planes": (type(None), list),
    "is3D": (type(None), bool),
    "remove_flat_cells": bool,
    "realtime_viewer": bool,
    "realtime_viewer_port": int,
    "realtime_viewer_max_cells": (type(None), int),
    "realtime_viewer_fixed_radius": (type(None), float),
    "mean_gene_counts_per_class": int,
    "mean_gene_counts_per_cell": int,
}


class Config(dict):
    """
    Simple dict-based configuration with validation.

    Inherits from dict for easy access and serialization, while adding
    validation and logging capabilities.

    Example:
        >>> cfg = Config({'max_iter': 500})
        >>> cfg['max_iter']
        500
        >>> cfg.set_runtime_attrs(coo_matrices)
    """

    def __init__(self, user_opts: Optional[Dict[str, Any]] = None):
        """
        Create config from defaults + user overrides.

        Args:
            user_opts: Optional dictionary of user configuration overrides.
                      Unknown keys will generate warnings but won't fail.

        Raises:
            TypeError: If any config parameter has incorrect type
        """
        # Start with defaults
        super().__init__(config.DEFAULT.copy())

        # Merge user options
        if user_opts:
            self._merge_user_opts(user_opts)

        # Validate types
        self._validate_types()

        # Setup logging
        self._setup_logging()

    def _merge_user_opts(self, opts: Dict[str, Any]) -> None:
        """
        Merge user options with defaults.

        Logs each override and warns about unrecognized keys.
        Unrecognized keys are still preserved (e.g., for callbacks).

        Args:
            opts: User configuration dictionary
        """
        valid_keys = set(config.DEFAULT.keys())

        for key, value in opts.items():
            if key in valid_keys:
                self[key] = value
                logger.info(f"Config override: {key} = {value}")
            else:
                # Preserve unknown keys (e.g., callbacks, plugins)
                self[key] = value
                logger.debug(f"Non-standard config key preserved: '{key}'")

    def _validate_types(self) -> None:
        """
        Validate all config parameters against type specifications.

        Raises:
            TypeError: If any parameter has incorrect type
        """
        for param_name, allowed_types in _TYPE_SPECS.items():
            # Skip if not in config (runtime params set later)
            if param_name not in self:
                continue

            value = self[param_name]

            # Normalize to tuple for isinstance()
            if not isinstance(allowed_types, tuple):
                allowed_types = (allowed_types,)

            # Check type
            if not isinstance(value, allowed_types):
                # Format error message
                type_names = " or ".join(
                    t.__name__ if hasattr(t, "__name__") else str(t)
                    for t in allowed_types
                )
                raise TypeError(
                    f"Config parameter '{param_name}' must be {type_names}, "
                    f"got {type(value).__name__}"
                )

    def _setup_logging(self) -> None:
        """Setup file handler for logging."""
        from pciSeq.src.core.utils.io_utils import log_file

        log_file(self)

    def set_runtime_attrs(self, coo) -> None:
        """
        Set configuration attributes that can only be determined at runtime.

        Args:
            coo: Either a single coo_matrix or list of coo_matrices

        Updates:
            - is3D: Whether data is 3D (multiple planes)
            - is_redis_running: Whether Redis server is available
            - exclude_planes: Normalized to empty list if None
        """
        from pciSeq.src.diagnostics.utils import check_redis_server

        self["is3D"] = self._detect_3d(coo)
        self["is_redis_running"] = check_redis_server()
        self["exclude_planes"] = self["exclude_planes"] or []

    def _detect_3d(self, coo) -> bool:
        """
        Detect if data is 3D based on coo structure.

        Args:
            coo: Either a coo_matrix or list of coo_matrices

        Returns:
            bool: True if 3D (multiple planes), False if 2D (single plane)

        Rules:
            - Single coo_matrix → 2D
            - List with one coo_matrix → 2D
            - List with multiple coo_matrices → 3D
        """
        from scipy.sparse import coo_matrix

        # Single matrix = 2D
        if isinstance(coo, coo_matrix):
            return False

        # List of matrices
        if isinstance(coo, list):
            if len(coo) == 1:
                return False  # Single plane = 2D
            elif len(coo) > 1:
                return True  # Multiple planes = 3D

        raise TypeError("coo must be coo_matrix or list of coo_matrices")
