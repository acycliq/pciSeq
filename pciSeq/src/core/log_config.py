import numpy as np
from pathlib import Path
import sys
import logging


def attach_to_log():
    """
    exists only for backwards compatibility.
    Replaced by logger_setup
    """
    logger_setup()


def logger_setup():
    """
    Taken from cellpose
    """
    cp_dir = Path.home().joinpath('.pciSeq')
    cp_dir.mkdir(exist_ok=True)
    log_file = cp_dir.joinpath('pciSeq.log')
    try:
        log_file.unlink()
    except:
        print('creating new log file')
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    return logger, log_file
