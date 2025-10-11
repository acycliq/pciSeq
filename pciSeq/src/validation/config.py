"""
Simplified configuration management for pciSeq.

This module provides a simpler dict-based Config class that replaces
the verbose 37-attribute dataclass ConfigManager.
"""
from typing import Dict, Any, Optional
from pciSeq import config
import logging

logger = logging.getLogger(__name__)


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
        """
        # Start with defaults
        super().__init__(config.DEFAULT.copy())

        # Merge user options
        if user_opts:
            self._merge_user_opts(user_opts)

        # Setup logging
        self._setup_logging()

    def _merge_user_opts(self, opts: Dict[str, Any]) -> None:
        """
        Merge user options with defaults.

        Logs each override and warns about unrecognized keys.

        Args:
            opts: User configuration dictionary
        """
        valid_keys = set(config.DEFAULT.keys())

        for key, value in opts.items():
            if key in valid_keys:
                self[key] = value
                logger.info(f'Config override: {key} = {value}')
            else:
                logger.warning(
                    f"Unrecognized config key: '{key}'. "
                    f"Valid keys: {', '.join(sorted(valid_keys))}"
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

        self['is3D'] = self._detect_3d(coo)
        self['is_redis_running'] = check_redis_server()
        self['exclude_planes'] = self['exclude_planes'] or []

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
                return True   # Multiple planes = 3D

        raise TypeError("coo must be coo_matrix or list of coo_matrices")