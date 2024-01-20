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

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    return logger
