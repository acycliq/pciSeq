from dataclasses import dataclass
from typing import List, Dict, Optional, Union
from numbers import Number
import logging
from pciSeq import config
from pciSeq.src.core.utils import log_file, check_redis_server
config_manager_logger = logging.getLogger(__name__)


@dataclass
class ConfigManager:
    exclude_genes: List[str]
    max_iter: int
    CellCallTolerance: float
    rGene: int
    Inefficiency: float
    InsideCellBonus: Union[bool, int]
    MisreadDensity: Union[float, Dict[str, float]]
    SpotReg: float
    nNeighbors: int
    rSpot: int
    save_data: bool
    output_path: List[str]
    launch_viewer: bool
    launch_diagnostics: bool
    is_redis_running: bool
    cell_radius: Optional[float]
    cell_type_prior: str
    mean_gene_counts_per_class: int
    mean_gene_counts_per_cell: int

    @classmethod
    def from_opts(cls, opts: Optional[Dict] = None) -> 'ConfigManager':
        """Create configuration from default values and optional overrides"""
        # Start with default configuration
        cfg_dict = config.DEFAULT.copy()

        log_file(cfg_dict)
        cfg_dict['is_redis_running'] = check_redis_server()

        # Override with user options if provided
        if opts is not None:
            default_items = set(cfg_dict.keys())
            user_items = set(opts.keys())
            if not user_items.issubset(default_items):
                raise ValueError(f'Options passed-in should be a dict with keys: {default_items}')

            for key, value in opts.items():
                if isinstance(value, (int, float, list, str, dict)) or \
                        (callable(getattr(value, '__call__', None)) and isinstance(value(1), Number)):
                    cfg_dict[key] = value
                else:
                    raise TypeError(
                        f"Invalid type for {key}. Only integers, floats, lists, strings, and dicts are allowed")
                config_manager_logger.info(f'{key} is set to {value}')

        # Create instance
        instance = cls(**cfg_dict)
        # instance._validate()
        return instance

    def to_dict(self) -> Dict:
        """Convert configuration back to dictionary format"""
        return {k: v for k, v in self.__dict__.items() if not k.startswith('_')}