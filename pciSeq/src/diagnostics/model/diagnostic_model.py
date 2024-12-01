"""
Model component of the MVC pattern for pciSeq diagnostics.
"""

import redis
import pandas as pd
import numpy as np
import logging
from typing import Optional, Dict, Any
from pciSeq.src.diagnostics.constants import DiagnosticKeys

model_logger = logging.getLogger(__name__)


class DiagnosticModel:
    def __init__(self, redis_host: str = 'localhost', redis_port: int = 6379):
        """Initialize Redis connection."""
        self.redis_client = redis.Redis(
            host=redis_host,
            port=redis_port,
            decode_responses=True
        )
        self._verify_connection()

    def publish_diagnostics(self, algorithm_model: Any, iteration: int, has_converged: bool) -> None:
        """Publish diagnostic data to Redis."""
        try:
            # Publish gene efficiency data
            gene_data = {
                'data': pd.DataFrame({
                    'gene_efficiency': algorithm_model.genes.eta_bar,
                    'gene': algorithm_model.genes.gene_panel
                }).to_json(),
                'metadata': {
                    'iteration': iteration,
                    'has_converged': has_converged,
                    'timestamp': pd.Timestamp.now().isoformat()
                }
            }
            self.redis_client.set(DiagnosticKeys.GENE_EFFICIENCY.value, str(gene_data))
            self.redis_client.publish(DiagnosticKeys.GENE_EFFICIENCY.value, 'update')

            # Publish cell type distribution data
            cell_type_indices = np.argmax(algorithm_model.cells.classProb[1:, :], axis=1)
            counts = np.bincount(
                cell_type_indices,
                minlength=len(algorithm_model.cellTypes.names)
            )

            cell_data = {
                'data': pd.DataFrame({
                    'class_name': algorithm_model.cellTypes.names,
                    'counts': counts
                }).to_json(),
                'metadata': {
                    'iteration': iteration,
                    'has_converged': has_converged,
                    'timestamp': pd.Timestamp.now().isoformat()
                }
            }
            self.redis_client.set(DiagnosticKeys.CELL_TYPE_POSTERIOR.value, str(cell_data))
            self.redis_client.publish(DiagnosticKeys.CELL_TYPE_POSTERIOR.value, 'update')
        except Exception as e:
            model_logger.error(f"Failed to publish diagnostics: {e}")

    def get_diagnostic_data(self, key: DiagnosticKeys) -> Optional[Dict]:
        """Retrieve data from Redis."""
        try:
            data = self.redis_client.get(key.value)
            return eval(data) if data else None
        except Exception as e:
            model_logger.error(f"Failed to retrieve data: {e}")
            return None

    def _verify_connection(self) -> None:
        """Verify Redis connection is working."""
        try:
            self.redis_client.ping()
        except redis.ConnectionError as e:
            model_logger.error(f"Redis connection failed: {e}")
            raise
