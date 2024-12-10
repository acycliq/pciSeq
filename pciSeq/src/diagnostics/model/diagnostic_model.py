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
            # Publish both gene efficiency and cell type data using the same pattern
            self._publish_data_to_redis(
                key=DiagnosticKeys.GENE_EFFICIENCY,
                data=self._prepare_gene_data(algorithm_model),
                metadata=self._create_metadata(iteration, has_converged)
            )

            self._publish_data_to_redis(
                key=DiagnosticKeys.CELL_TYPE_POSTERIOR,
                data=self._prepare_cell_type_data(algorithm_model),
                metadata=self._create_metadata(iteration, has_converged)
            )
        except Exception as e:
            model_logger.error(f"Failed to publish diagnostics: {e}")

    def _publish_data_to_redis(self, key: DiagnosticKeys, data: pd.DataFrame, metadata: dict) -> None:
        """Helper method to publish data to Redis with a consistent pattern."""
        formatted_data = {
            'data': data.to_json(),
            'metadata': metadata
        }
        # Store the data in Redis
        self.redis_client.set(key.value, str(formatted_data))
        # Notify subscribers that new data are available
        self.redis_client.publish(key.value, 'update')

    @staticmethod
    def _prepare_gene_data(algorithm_model: Any) -> pd.DataFrame:
        """Prepare gene efficiency data."""
        return pd.DataFrame({
            'gene_efficiency': algorithm_model.genes.eta_bar,
            'gene': algorithm_model.genes.gene_panel
        })

    @staticmethod
    def _prepare_cell_type_data(algorithm_model: Any) -> pd.DataFrame:
        """Prepare cell type distribution data."""
        cell_type_indices = np.argmax(algorithm_model.cells.classProb[1:, :], axis=1)
        counts = np.bincount(
            cell_type_indices,
            minlength=len(algorithm_model.cellTypes.names)
        )
        return pd.DataFrame({
            'class_name': algorithm_model.cellTypes.names,
            'counts': counts
        })

    @staticmethod
    def _create_metadata(iteration: int, has_converged: bool) -> dict:
        """Create consistent metadata structure."""
        return {
            'iteration': iteration,
            'has_converged': has_converged,
            'timestamp': pd.Timestamp.now().isoformat()
        }

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

    def flush_db(self):
        """Clear all keys from the Redis database."""
        self.redis_client.flushdb()
