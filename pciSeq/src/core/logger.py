import sys
import logging


def attach_to_log():
    """
    exists only for backwards compatibility.
    Replaced by logger_setup
    """
    setup_logger()


def setup_logger():
    # Remove all handlers associated with the root logger object.
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    return logger
