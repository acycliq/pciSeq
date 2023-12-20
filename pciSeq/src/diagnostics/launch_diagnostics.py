import os
import sys
import subprocess
from pciSeq.src.core.log_config import logger, attach_to_log
import pciSeq.src.diagnostics.utils as utils


def launch_dashboard():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')
    # redirect to dev/null to bypass the welcome message from streamlit
    # devnull = '> /dev/null' if sys.platform in ['linux', 'darwin'] else '> NUL'
    code_1, code_2 = utils.validate_redis()
    if code_1 == 0 and code_2 == 0:
        p = subprocess.Popen(["streamlit", "run", filename, os.devnull])
        # TODO: you need to kill the process on exit
        # logger.info('Starting process with pid: %d to run the diagnostics' % p.pid)
    else:
        logger.info("Skipping diagnostics, cannot run them. Either redis not installed or not running.")


if __name__ == "__main__":
    attach_to_log()
    launch_dashboard()
