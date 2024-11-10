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
    InsideCellBonus: Union[bool, int, float]
    MisreadDensity: Union[float, Dict[str, float]]
    SpotReg: float
    nNeighbors: int
    rSpot: int
    save_data: bool
    output_path: str
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
        if opts is None:
            opts = config.DEFAULT.copy()

        # Start with default configuration
        cfg_dict = config.DEFAULT.copy()

        # # Override with user options if provided
        for key in opts:
            if key in cfg_dict:
                cfg_dict[key] = opts[key]
                config_manager_logger.info(f'{key} is set to {opts[key]}')

        log_file(cfg_dict)
        cfg_dict['is_redis_running'] = check_redis_server()

        # Create instance
        instance = cls(**cfg_dict)
        # instance._validate()
        return instance

    def to_dict(self) -> Dict:
        """Convert configuration back to dictionary format"""
        return {k: v for k, v in self.__dict__.items() if not k.startswith('_')}