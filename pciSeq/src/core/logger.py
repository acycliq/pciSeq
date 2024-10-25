import sys
import logging
import colorlog


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
        level = logging.INFO

    # Create a logger specific to pciSeq
    logger = logging.getLogger('pciSeq')

    if not logger.handlers:
        logger.setLevel(level)
        logger.propagate = True  # Allow propagation to root logger

        # Create color formatter
        color_formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            log_colors={
                'DEBUG': 'cyan',
                'INFO': 'green',
                'WARNING': 'yellow',
                'ERROR': 'red',
                'CRITICAL': 'red,bg_white',
            },
            reset=True,
            style='%'
        )

        # Console handler
        console_handler = colorlog.StreamHandler(sys.stdout)
        console_handler.setFormatter(color_formatter)
        logger.addHandler(console_handler)

        # File handler (optional, without colors)
        if log_file:
            file_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
            logger.info('Writing to %s' % log_file)

    return logger


def get_logger(name):
    """
    Get a logger that is a child of the main pciSeq logger.

    :param name: Name of the module requesting the logger
    :return: Logger instance
    """
    logger = logging.getLogger(f'pciSeq.{name}')
    return logger