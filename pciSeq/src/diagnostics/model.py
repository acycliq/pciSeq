import pciSeq.src.diagnostics.utils as utils
import pandas as pd
import pickle


class DiagnosticsModel:
    def __init__(self):
        self.redis_client = utils.redis_db(flush=False).redis_client

    def is_redis_available(self):
        return utils.validate_redis() == (0, 0)

    def get_gene_efficiency(self):
        data = self.redis_client.get("gene_efficiency")
        return pickle.loads(data) if data else None

    def get_cell_type_posterior(self):
        data = self.redis_client.get("cell_type_posterior")
        return pickle.loads(data) if data else None

    def subscribe_to_channels(self):
        p = self.redis_client.pubsub()
        subscribe_to = ['gene_efficiency', 'cell_type_posterior']
        p.psubscribe(*subscribe_to)
        return p
