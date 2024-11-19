import pandas as pd
import numpy as np
from .constants import DiagnosticKeys
import time

# This module is part of the diagnostics package in the pciSeq project.
# It implements the RedisPublisher class, which is responsible for
# publishing diagnostic information to a Redis database.

# In the context of the Model-View-Controller (MVC) pattern:
# - This class acts primarily as part of the Controller.
# - It updates the Model (Redis database) with the latest algorithm state.
# - The View would be a separate component that reads from Redis to display data.
#
# Note: The actual View component is not present in this file.
# It would typically be a separate module or application that
# subscribes to Redis and presents the data to the user.
# The View component is primarily implemented in the diagnostics.py
# file within the pciSeq project.


class RedisPublisher:
    def __init__(self, redis_db):
        self.redis_db = redis_db

    def publish_diagnostics(self, model, iter_num, has_converged):
        """
        Publish diagnostic information to Redis.
        This method acts as part of the Controller in the MVC pattern,
        updating the Model (Redis database) with the latest algorithm state.
        """
        self.publish_gene_efficiency_to_redis(model, iter_num, has_converged)
        self.publish_cell_type_distribution_to_redis(model, iter_num, has_converged)

    def publish_gene_efficiency_to_redis(self, model, iter_num, has_converged):
        """
        Publish gene efficiency information to Redis.
        """
        gene_efficiency_data = pd.DataFrame({
            'gene_efficiency': model.genes.eta_bar,
            'gene': model.genes.gene_panel
        })
        self._publish_to_redis(gene_efficiency_data, DiagnosticKeys.GENE_EFFICIENCY.value, iter_num, has_converged)

    def publish_cell_type_distribution_to_redis(self, model, iter_num, has_converged):
        """
        Publish cell type distribution information to Redis.
        """
        cell_type_counts = self._get_cell_type_counts(model)
        cell_type_data = pd.DataFrame({
            'class_name': model.cellTypes.names,
            'counts': cell_type_counts
        })
        self._publish_to_redis(cell_type_data, DiagnosticKeys.CELL_TYPE_POSTERIOR.value, iter_num, has_converged)

    def _get_cell_type_counts(self, model):
        """
        Get the count of cells assigned to each cell type.
        """
        cell_type_indices = np.argmax(model.cells.classProb[1:, :], axis=1)
        return np.bincount(cell_type_indices, minlength=len(model.cellTypes.names))

    def _publish_to_redis(self, data, key, iter_num, has_converged):
        """
        Publish data to Redis with metadata.
        This method is the final step in updating the Model (Redis) with new data.
        """
        if self.redis_db:
            self.redis_db.publish(data, key,
                                  iteration=iter_num,
                                  has_converged=has_converged,
                                  unix_time=time.time())


