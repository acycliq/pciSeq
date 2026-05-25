import subprocess
import sys
import os
from pciSeq._version import __version__


def _resolve_git_info():
    """Return (commit, branch). Try live git first, fall back to values
    baked at build time, fall back to 'unknown'."""
    pkg_dir = os.path.dirname(__file__)

    def _git(args):
        try:
            return subprocess.check_output(
                ['git'] + args,
                cwd=pkg_dir,
                stderr=subprocess.DEVNULL,
                text=True,
            ).strip()
        except Exception:
            return ''

    commit = _git(['rev-parse', '--short', 'HEAD'])
    branch = _git(['rev-parse', '--abbrev-ref', 'HEAD'])
    if commit and branch:
        return commit, branch

    try:
        from pciSeq._build_info import __commit__ as c, __branch__ as b
        return c or 'unknown', b or 'unknown'
    except ImportError:
        return 'unknown', 'unknown'


def _resolve_build_date():
    """Date setup.py was last run (write time of _build_info.py).
    'unknown' if the package was imported from a source checkout that
    has never been built."""
    try:
        from pciSeq._build_info import __build_date__
        return __build_date__ or 'unknown'
    except ImportError:
        return 'unknown'


__commit__, __branch__ = _resolve_git_info()
__build_date__ = _resolve_build_date()


from pciSeq.app import fit
from pciSeq.app import cell_type
from pciSeq.src.preprocess.main import stage_data
from pciSeq.src.core.logger import attach_to_log, setup_logger
# from pciSeq.src.core.analysis import CellExplorer
import logging

logger = logging.getLogger(__name__)


def confirm_prompt(question):
    reply = None
    while reply not in ("", "y", "n"):
        reply = input(f"{question} (y/n): ").lower()
    return reply in ("", "y")


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
    from pciSeq.src.viewer.stage_image import tile_maker, stage_image
else:
    def tile_maker():
        logger.warning('>>>> tile_maker() isnt available because libvips is not installed. Please see '
                            'https://www.libvips.org/install.html <<<<')
        logger.warning('>>>> If you are on Linux you can install it by calling: sudo apt install libvips <<<<')

    def stage_image(*args, **kwargs):
        logger.warning('>>>> stage_image() isnt available because libvips is not installed. Please see '
                            'https://www.libvips.org/install.html <<<<')
        logger.warning('>>>> If you are on Linux you can install it by calling: sudo apt install libvips <<<<')
