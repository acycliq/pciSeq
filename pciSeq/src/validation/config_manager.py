from dataclasses import dataclass
from typing import List, Dict, Optional, Union
import logging
from scipy.sparse import coo_matrix
from pciSeq import config
from pciSeq.src.core.utils.io_utils import log_file
from pciSeq.src.diagnostics.utils import check_redis_server
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
    cell_centroid_prior: Union[float, int, Dict[str, float]]
    cell_cov_prior: Union[float, int, Dict[str, float]]
    SpotReg: float
    nNeighbors: int
    rSpot: Union[int, float]
    save_data: bool
    output_path: str
    launch_viewer: Union[bool, str]
    launch_diagnostics: bool
    is_redis_running: bool
    cell_radius: Optional[float]
    cell_type_prior: str
    voxel_size: list
    exclude_planes: list
    is3D: Union[None, bool]
    remove_flat_cells: bool
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
            else:
                config_manager_logger.warning(f"Unrecognized configuration option: '{key}'! "
                                              f"Valid options are: {', '.join(sorted(cfg_dict.keys()))}")

        log_file(cfg_dict)

        # Create instance
        instance = cls(**cfg_dict)
        # instance._validate()
        return instance

    def to_dict(self) -> Dict:
        """Convert configuration back to dictionary format"""
        return {k: v for k, v in self.__dict__.items() if not k.startswith('_')}

    def set_runtime_attributes(self, coo):
        """Set configuration attributes that can only be determined at runtime.

        Args:
            coo: Coordinate data

        Returns:
            config: Updated configuration with runtime attributes
        """
        self.is3D = self.check_is3D(coo)
        self.is_redis_running = check_redis_server()

        # if exclude_planes is None set it to []
        self.exclude_planes = self.exclude_planes or []

    def check_is3D(self, input_data: Union[coo_matrix, List[coo_matrix]]) -> bool:
        """
        Determine if the input is a single COO matrix (is3D=False) or multiple COO matrices (is3D=True).

        Args:
            input_data: Either a single coo_matrix or a list of coo_matrices

        Returns:
            bool: False if input is a single coo_matrix, True if it's a list with multiple coo_matrices
        """
        if isinstance(input_data, coo_matrix):
            return False
        elif isinstance(input_data, list):
            if len(input_data) == 1 and isinstance(input_data[0], coo_matrix):
                return False
            elif all(isinstance(mat, coo_matrix) for mat in input_data):
                return True
            else:
                raise ValueError("If input is a list, all elements must be coo_matrix")
        else:
            raise TypeError("Input must be either a coo_matrix or a list of coo_matrices")
