from enum import Enum, auto


class DiagnosticKeys(Enum):
    """
    This enumeration serves as the centralised source for all keys used in
    Redis publish/subscribe operations throughout the diagnostics package:
        - RedisPublisher uses these keys to publish diagnostic data
        - DiagnosticsModel uses these keys to subscribe to and retrieve data

    Current keys:
        - GENE_EFFICIENCY: Used for gene efficiency metrics
        - CELL_TYPE_POSTERIOR: Used for cell type distribution data
    """
    GENE_EFFICIENCY = "gene_efficiency"
    CELL_TYPE_POSTERIOR = "cell_type_posterior"
