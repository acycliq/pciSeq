"""Writing a finished run to disk.

Split by format: tsv_export for the plain text tables, arrow_export for the
feather files the viewer streams, diagnostics_db for the sqlite database it
queries, spatialdata_export (over in utils) for the zarr store. export.py is
the one function that calls them all.

Everything a caller outside this package needs is re-exported here, so
`from pciSeq.src.core.io import write_data` keeps working no matter how the
modules underneath get shuffled around later.
"""

from .export import write_data
from .paths import get_out_dir, get_pciSeq_install_dir, log_file
from .remote import download_url_to_file, load_from_url
from .provenance import run_metadata, serialise
from .tsv_export import write_tsv, read_tsv
from .arrow_export import write_arrow, geneData_to_arrow, cellData_to_arrow, boundaries_to_arrow
from .diagnostics_db import export_diagnostics, shrink_diagnostics_db

__all__ = [
    'write_data',
    'get_out_dir', 'get_pciSeq_install_dir', 'log_file',
    'download_url_to_file', 'load_from_url',
    'run_metadata', 'serialise',
    'write_tsv', 'read_tsv',
    'write_arrow', 'geneData_to_arrow', 'cellData_to_arrow', 'boundaries_to_arrow',
    'export_diagnostics',
    'shrink_diagnostics_db',
]
