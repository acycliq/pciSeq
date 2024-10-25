import subprocess
import sys
import os
from pciSeq._version import __version__
from pciSeq.app import fit
from pciSeq.app import cell_type
from pciSeq.src.preprocess.spot_labels import stage_data
import pciSeq.src.core.utils as utils
from pciSeq.src.core.logger import attach_to_log, setup_logger
from pciSeq.src.core.logger import get_logger
import logging

# Add a NullHandler to the pciSeq logger
logging.getLogger('pciSeq').addHandler(logging.NullHandler())

# the following creates a logger named pciSeq.pciSeq (which is a child of the pciSeq logger)
logger = get_logger(__name__)

# Some notes on logging:
# 1. Logger Propagation: By default, loggers propagate their messages up the hierarchy to their
#    parent loggers, ultimately reaching the root logger.
# 2. NullHandler behavior: The NullHandler we added to the 'pciSeq' logger doesn't actually
#    process or output any log messages. It simply exists to prevent "No handlers could be found"
#    warnings.
#    Here's how it works:
#    * We add a NullHandler to the 'pciSeq' logger:
#    * This 'pciSeq' logger now has a handler, so it won't raise
#      "No handlers could be found" warnings.
# 3. However, because the NullHandler doesn't actually do anything with the log messages,
#    and because logger propagation is still enabled (it's enabled by default), any log
#    messages sent to the 'pciSeq' logger will still propagate up to its parent loggers,
#    eventually reaching the root logger.
# 4. If an application using pciSeq has configured the root logger (or any parent logger of
#    'pciSeq'), those configurations will be used to handle the log messages that have
#    propagated up from the 'pciSeq' logger.
#
# So, we're "still allowing the root logger to handle the actual logging output" by:
# 1. Not disabling propagation on the 'pciSeq' logger.
# 2. Not adding any handlers to the 'pciSeq' logger that would actually process
#    or output log messages (the NullHandler doesn't do this).
# This setup ensures that pciSeq doesn't interfere with the logging configuration
#   of applications that use it, while also preventing warnings about missing handlers.


def confirm_prompt(question):
    reply = None
    while reply not in ("", "y", "n"):
        reply = input(f"{question} (y/n): ").lower()
    return (reply in ("", "y"))


def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


def install_libvips():
    subprocess.check_call("apt-get update", shell=True)
    subprocess.check_call("apt-get install", shell=True)
    subprocess.check_call(['apt-get', 'install', '-y', 'libvips'],
               stdout=open(os.devnull, 'wb'), stderr=subprocess.STDOUT)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyvips"])

#
# def check_libvips(logger):
#     confirm = confirm_prompt('Install libvips?')
#     if confirm:
#       install_libvips()
#     else:
#       print('>>>> libvips not installed')
#     return confirm


def check_libvips():
    try:
        import pyvips
        status = True
    except OSError:
        status = False
    except Exception as err:
        raise
    return status


if check_libvips():
    from pciSeq.src.viewer.stage_image import tile_maker
else:
    def tile_maker():
        logger.warning('>>>> tile_maker() because libvips is not installed. Please see https://www.libvips.org/install.html <<<<')
        logger.warning('>>>> If you are on Linux you can install it by calling: sudo apt install libvips <<<<')



