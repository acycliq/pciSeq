import logging
import sys


def setup_logger(name):
    """
    Set up a logger that writes to standard output
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # Check if logger already has handlers to avoid duplicates
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger
