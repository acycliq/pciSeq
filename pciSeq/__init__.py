import subprocess
import sys
import os
from pciSeq._version import __version__
from pciSeq.app import fit
from pciSeq.app import cell_type
from pciSeq.src.preprocess.spot_labels import stage_data
import pciSeq.src.core.utils as utils
from pciSeq.src.core.logger import attach_to_log, setup_logger
from pciSeq.src.core.analysis import CellExplorer
import logging

init_logger = logging.getLogger(__name__)


def confirm_prompt(question):
    reply = None
    while reply not in ("", "y", "n"):
        reply = input(f"{question} (y/n): ").lower()
    return (reply in ("", "y"))


def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


def install_libvips():
    try:
        subprocess.check_call("apt-get update", shell=True)
        subprocess.check_call("apt-get install -y libvips", shell=True)  # Combined command with -y flag
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyvips"])
        return True
    except subprocess.CalledProcessError as e:
        init_logger.error(f"Failed to install libvips: {str(e)}")
        return False

def check_libvips():
    try:
        import pyvips
        return True
    except ImportError:  # More specific exception
        init_logger.warning('libvips not found. Please install it manually or run with sudo privileges.')
        return False
    except Exception as err:
        init_logger.error(f"Unexpected error checking libvips: {str(err)}")
        raise


if check_libvips():
    from pciSeq.src.viewer.stage_image import tile_maker
else:
    def tile_maker():
        init_logger.warning('>>>> tile_maker() because libvips is not installed. Please see https://www.libvips.org/install.html <<<<')
        init_logger.warning('>>>> If you are on Linux you can install it by calling: sudo apt install libvips <<<<')



