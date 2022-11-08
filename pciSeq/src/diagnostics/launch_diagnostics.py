import os
import sys
import subprocess
from pciSeq.src.cell_call.log_config import logger


def launch_dashboard():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')
    exe = sys.executable

    p = subprocess.Popen([exe, "-m" "streamlit", "run", filename])
    # logger.info('Starting process with pid: %d to run the diagnostics' % p.pid)


if __name__ == "__main__":
    launch_dashboard()
