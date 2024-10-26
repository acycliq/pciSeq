import sys
import colorlog
import logging


def attach_to_log():
    """
    exists only for backwards compatibility.
    Replaced by logger_setup
    """
    setup_logger()


def setup_logger(level=None):

    if level is None:
        level = logging.INFO

    # Remove all handlers associated with the root logger object.
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Create color formatter
    color_formatter = colorlog.ColoredFormatter(
        "%(log_color)s%(asctime)s - %(levelname)s - %(message)s",
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

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            console_handler
        ]
    )
    logger = logging.getLogger('pciSeq')

    return logger
