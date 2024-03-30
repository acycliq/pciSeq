import os
from pathlib import Path
import tomlkit
import subprocess
import pciSeq.src.diagnostics.utils as utils
import logging

launch_diagnostics_logger = logging.getLogger(__name__)


def launch_dashboard():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')
    code_1, code_2 = utils.validate_redis()
    if code_1 == 0 and code_2 == 0:
        make_credentials()
        p = subprocess.Popen(["streamlit", "run", filename, " --server.headless true"])
        # TODO: you need to kill the process on exit
        # launch_diagnostics_logger.info('Starting process with pid: %d to run the diagnostics' % p.pid)
    else:
        launch_diagnostics_logger.info("Skipping diagnostics, cannot run them. Either redis not installed or not running.")


def make_credentials():
    """
    creates a credentials.toml file in the .streamlit folder under the user's home dir and
    in the [general] section it creates an email key with value an empty string.
    If an email key/value pair exists already, the file remains intact.
    The only purpose of this is to skip the annoying streamlit welcome message so that
    the diagnostics will be launched without any user interaction. (otherwise the welcome msg
    might pause the program flow)
    """
    credentials = os.path.join(Path.home(), '.streamlit', 'credentials.toml')
    Path(credentials).parent.mkdir(parents=True, exist_ok=True)
    mode = 'rt' if os.path.exists(credentials) else 'w+'
    with open(credentials, mode) as fp:
        doc = tomlkit.load(fp)
        append_email(doc, credentials)


def append_email(doc, fname):
    if doc.get('general') is None:
        general = tomlkit.table()
        doc.add("general", general)
    if doc.get('general').get('email') is None:
        doc.get('general').add("email", "")
    with open(fname, 'w') as outfile:
        tomlkit.dump(doc, outfile)