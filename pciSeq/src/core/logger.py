import sys
import logging


def attach_to_log():
    """
    exists only for backwards compatibility.
    Replaced by logger_setup
    """
    setup_logger()


def setup_logger(level=None, log_file=None):
    """
    Set up the pciSeq logger. Should be called once

    :param level: The logging level (default: logging.INFO)
    :param log_file: Path to a log file (optional)
    :return: Configured logger
    """

    if level is None:
        level=logging.INFO

    # Create a logger specific to pciSeq
    logger = logging.getLogger('pciSeq')

    # Only configure if the logger doesn't already have handlers
    if not logger.handlers:
        logger.setLevel(level)

        # Create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # File handler (optional)
        if log_file:
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            logger.info('Writing to %s' % log_file)

        # Prevent propagation to root logger
        logger.propagate = False

    return logger


def get_logger(name):
    """
    Get a logger that is a child of the main pciSeq logger.

    :param name: Name of the module requesting the logger
    :return: Logger instance
    """
    return logging.getLogger(f'pciSeq.{name}')