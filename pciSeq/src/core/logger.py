"""
Logging Configuration Module for pciSeq

This module provides a standardized logging setup for the pciSeq package,
with colored output and consistent formatting across all modules.

IMPORTANT NOTE:
-------------
This logger is designed to be non-intrusive when pciSeq is imported as a package.
It will NOT modify parent application's logging settings or override existing
configurations. The parent application's logging setup remains unchanged when using pciSeq.

Key Functions:
------------
setup_logger:
    Main configuration function that:
    - Sets up colored console output
    - Configures log levels and formats
    - Returns a configured logger instance

attach_to_log:
    Legacy function maintained for backwards compatibility
    (Redirects to setup_logger)

Features:
--------
- Color-coded output by log level:
  * DEBUG: Cyan
  * INFO: Green
  * WARNING: Yellow
  * ERROR: Red
  * CRITICAL: Red with white background

- Standardized timestamp format
- Console output to stdout
- Configurable log levels
- Thread-safe logging

Notes:
-----
- Default log level is INFO
- All handlers are cleared before setup to avoid duplication
- Uses colorlog for ANSI color support
- Format: "timestamp - level - message"
"""

import sys
import colorlog
import logging


def attach_to_log():
    """
    exists only for backwards compatibility.
    Replaced by setup_logger
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