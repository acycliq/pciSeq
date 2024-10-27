import pciSeq.src.diagnostics.utils as utils
import pickle

"""
This module contains the DiagnosticsModel class, which serves as the Model component
in the Model-View-Controller (MVC) pattern for the diagnostics functionality.

The DiagnosticsModel is responsible for managing the data and business logic
related to diagnostics in the pciSeq application. It interacts with the Redis
database to store and retrieve diagnostic information, and provides methods
for accessing this data.

Key responsibilities of the DiagnosticsModel include:
1. Establishing and maintaining a connection to the Redis database.
2. Checking the availability of the Redis server.
3. Retrieving specific diagnostic data (e.g., gene efficiency, cell type posterior).
4. Setting up subscriptions to Redis channels for real-time updates.

This class acts as an intermediary between the Redis database and the rest of the
application, encapsulating the data access logic and providing a clean interface
for the Controller to interact with the Model.
"""


def is_redis_available():
    """Check if Redis server is available and running."""
    return utils.validate_redis() == (0, 0)


class DiagnosticsModel:
    """
    Manages diagnostic data interactions with Redis for pciSeq.

    Attributes:
        redis_client (redis.Redis): Redis client instance for database operations.
    """
    def __init__(self):
        self.redis_client = utils.RedisDB(flush=False).redis_client

    def get_gene_efficiency(self):
        """Retrieve gene efficiency data from Redis."""
        data = self.redis_client.get("gene_efficiency")
        return pickle.loads(data) if data else None

    def get_cell_type_posterior(self):
        """Retrieve cell type posterior data from Redis."""
        data = self.redis_client.get("cell_type_posterior")
        return pickle.loads(data) if data else None

    def subscribe_to_channels(self):
        """Set up subscriptions to relevant Redis channels."""
        p = self.redis_client.pubsub()
        subscribe_to = ['gene_efficiency', 'cell_type_posterior']
        p.psubscribe(*subscribe_to)
        return p
