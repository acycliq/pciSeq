"""What produced a run: the stamp that goes on it, and the pickled model."""

import os
import pickle
import logging
from typing import Dict, Any

from .diagnostics_db import export_diagnostics

logger = logging.getLogger(__name__)



def run_metadata() -> Dict:
    """The stamp that goes on every run, so a result can be traced back to the code
    and the environment that made it.

    VarBayes calls this once, when the model is built, and keeps the result as
    `metadata`. The pickle, the diagnostics db and the SpatialData store all carry
    that same dict, so there is one stamp and one place to add a field.

    Everything in it is a string, or a dict of strings, because the db and the store
    both write it out as json. There is no hostname or user name in here on purpose:
    these files get passed on to other people, and the machine a run happened on says
    nothing about how to reproduce it.
    """
    import sys
    import platform
    import datetime
    from importlib import metadata as importlib_metadata

    # imported here and not at the top: pciSeq/__init__.py is still half way through
    # its own imports when this module gets loaded
    from pciSeq import __version__, __branch__, __commit__, __commit_date__

    # the libraries the numbers depend on. numba is in because spots_to_cell runs
    # through a numba kernel.
    package_versions = {}
    for pkg in ('numpy', 'scipy', 'pandas', 'numba'):
        try:
            package_versions[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            package_versions[pkg] = 'unknown'

    return {
        'version': __version__,
        'branch': __branch__,
        'commit': __commit__,
        'commit_date': __commit_date__,
        'created_at': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'python_version': sys.version.split()[0],
        'os': f'{platform.system()} {platform.release()}',
        'package_versions': package_versions,
    }


def serialise(varBayes: Any, debug_dir: str) -> None:
    """Write diagnostics.db, then pickle variable Bayes object to debug directory.

    Args:
        varBayes: Object to serialize
        debug_dir: Directory to save pickle file
    """
    # the db goes first. The pickle is the big slow one, so if the save dies half
    # way the viewer still has its db.
    # the viewer reads everything under viewer_data, the db included
    viewer_dir = os.path.join(os.path.dirname(debug_dir), 'viewer_data')
    export_diagnostics(varBayes, viewer_dir)

    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)
    pickle_dst = os.path.join(debug_dir, 'pciSeq.pickle')
    with open(pickle_dst, 'wb') as outf:
        pickle.dump(varBayes, outf)

    pickle_mb = os.path.getsize(pickle_dst) / (1024 * 1024)
    logger.info('Saved at %s (%.1f MB)', pickle_dst, pickle_mb)
