import os
from pathlib import Path
import tomlkit
import subprocess
import pciSeq.src.diagnostics.utils as utils
import logging
from pciSeq.src.diagnostics.model import DiagnosticsModel

launch_logger = logging.getLogger(__name__)


class DiagnosticsController:
    def __init__(self):
        self.model = DiagnosticsModel()

    def launch_dashboard(self):
        """
        Launches a diagnostic dashboard using Streamlit if Redis is available.
        This function acts as a controller in the MVC pattern, coordinating the
        data access (model) and the dashboard display (view).
        """
        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, 'diagnostics.py')
        if self.model.is_redis_available():
            self.make_credentials()
            try:
                # The Streamlit script (diagnostics.py) should handle data presentation
                # and act as the view component
                process = subprocess.Popen(["streamlit", "run", filename, " --server.headless true"])
                launch_logger.info(f'Starting diagnostics dashboard with PID: {process.pid}')
                return process
            except subprocess.SubprocessError as e:
                launch_logger.error(f"Failed to start Streamlit: {e}")
        else:
            launch_logger.warning("Skipping diagnostics: Redis is not available.")
        return None

    def shutdown_dashboard(self, process):
        if process:
            process.terminate()
            process.wait()
            launch_logger.info(f'Terminated diagnostics dashboard with PID: {process.pid}')

    @staticmethod
    def make_credentials():
        """
        Creates a credentials.toml file in the .streamlit folder under the user's home dir and
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
            DiagnosticsController.append_email(doc, credentials)

    @staticmethod
    def append_email(doc, fname):
        if doc.get('general') is None:
            general = tomlkit.table()
            doc.add("general", general)
        if doc.get('general').get('email') is None:
            doc.get('general').add("email", "")
        with open(fname, 'w') as outfile:
            tomlkit.dump(doc, outfile)


def launch_dashboard():
    controller = DiagnosticsController()
    dashboard_process = controller.launch_dashboard()
    # ... (other main logic)
    # When done:
    # controller.shutdown_dashboard(dashboard_process)


# Usage
if __name__ == "__main__":
    launch_dashboard()
