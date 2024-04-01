import subprocess
import sys
import os
from pciSeq._version import __version__
from pciSeq.app import fit
from pciSeq.app import cell_type
from pciSeq.src.preprocess.spot_labels import stage_data
import pciSeq.src.core.utils as utils
from pciSeq.src.core.logger import attach_to_log, setup_logger
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
        init_logger.warning('>>>> tile_maker() because libvips is not installed. Please see https://www.libvips.org/install.html <<<<')
        init_logger.warning('>>>> If you are on Linux you can install it by calling: sudo apt install libvips <<<<')



