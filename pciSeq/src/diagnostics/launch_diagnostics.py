import os
import subprocess
import pciSeq.src.diagnostics.utils as utils
import logging

launch_diagnostics_logger = logging.getLogger(__name__)


def launch_dashboard():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')
    code_1, code_2 = utils.validate_redis()
    if code_1 == 0 and code_2 == 0:
        p = subprocess.Popen(["streamlit", "run", filename, " --server.headless true"])
        # TODO: you need to kill the process on exit
        # logger.info('Starting process with pid: %d to run the diagnostics' % p.pid)
    else:
        launch_diagnostics_logger.info("Skipping diagnostics, cannot run them. Either redis not installed or not running.")

