"""
Controller component of the MVC pattern for pciSeq diagnostics.

This module implements the Controller layer, responsible for:
1. Coordinating between Model and View components
2. Managing dashboard lifecycle (start/stop)
3. Handling data flow from algorithm to diagnostics
4. System initialization and shutdown

The Controller acts as the central coordinator, managing the flow of data
and control between the Model and View components.
"""

import os
import subprocess
import logging
import signal
from pathlib import Path
from typing import Optional
import tomlkit

from ..model.diagnostic_model import DiagnosticModel
from ..constants import DiagnosticKeys

controller_logger = logging.getLogger(__name__)


class DiagnosticController:
    """
    Main controller class coordinating the diagnostics system.

    This class manages the interaction between the Model (diagnostic_model.py)
    and View (dashboard.py) components, handling system lifecycle and data flow.
    """

    def __init__(self):
        """Initialize controller with required components."""
        self.model = DiagnosticModel()
        self.dashboard_process: Optional[subprocess.Popen] = None
        self._setup_signal_handlers()

    def _setup_signal_handlers(self) -> None:
        """Setup handlers for clean shutdown on system signals."""
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def _signal_handler(self, signum: int, frame) -> None:
        """
        Handle system signals for clean shutdown.

        Args:
            signum: Signal number
            frame: Current stack frame
        """
        controller_logger.info(f"Received signal {signum}, initiating shutdown...")
        self.shutdown()

    def initialize(self) -> bool:
        """
        Initialize the diagnostics system.

        Returns:
            bool: True if initialization successful, False otherwise
        """
        try:
            self._setup_streamlit_credentials()
            controller_logger.info("Diagnostics system initialized successfully")
            return True
        except Exception as e:
            controller_logger.error(f"Failed to initialize diagnostics: {e}")
            return False

    def launch_dashboard(self) -> bool:
        """
        Launch the Streamlit dashboard (View component).

        Returns:
            bool: True if dashboard launched successfully, False otherwise
        """
        if not self.initialize():
            return False

        try:
            # Get path to dashboard script
            dirname = os.path.dirname(os.path.dirname(__file__))
            dashboard_path = os.path.join(dirname, 'view', 'dashboard.py')

            # Launch dashboard process
            self.dashboard_process = subprocess.Popen([
                "streamlit", "run", dashboard_path, " --server.headless true"
            ])

            controller_logger.info(f'Started dashboard with PID: {self.dashboard_process.pid}')
            return True
        except subprocess.SubprocessError as e:
            controller_logger.error(f"Failed to start dashboard: {e}")
            return False

    def update_diagnostics(self, algorithm_model, iteration: int, has_converged: bool) -> None:
        """
        Update diagnostic data through the Model.

        Args:
            algorithm_model: Current state of the algorithm
            iteration: Current iteration number
            has_converged: Whether the algorithm has converged
        """
        try:
            self.model.publish_diagnostics(
                algorithm_model=algorithm_model,
                iteration=iteration,
                has_converged=has_converged
            )
            controller_logger.debug(f"Updated diagnostics for iteration {iteration}")
        except Exception as e:
            controller_logger.error(f"Failed to update diagnostics: {e}")

    def shutdown(self) -> None:
        """Clean shutdown of all components."""
        if self.dashboard_process:
            try:
                self.dashboard_process.terminate()
                self.dashboard_process.wait(timeout=5)
                controller_logger.info(f'Terminated dashboard with PID: {self.dashboard_process.pid}')
            except subprocess.TimeoutExpired:
                self.dashboard_process.kill()
                controller_logger.warning(f'Forced dashboard termination with PID: {self.dashboard_process.pid}')
            finally:
                self.dashboard_process = None

    def _setup_streamlit_credentials(self) -> None:
        """
        Setup Streamlit credentials to bypass welcome message.

        Creates or updates .streamlit/credentials.toml to skip the welcome screen.
        """
        credentials_path = os.path.join(Path.home(), '.streamlit', 'credentials.toml')
        Path(credentials_path).parent.mkdir(parents=True, exist_ok=True)

        try:
            mode = 'rt' if os.path.exists(credentials_path) else 'w+'
            with open(credentials_path, mode) as fp:
                doc = tomlkit.load(fp) if mode == 'rt' else tomlkit.document()

                if 'general' not in doc:
                    doc['general'] = tomlkit.table()
                if 'email' not in doc['general']:
                    doc['general']['email'] = ""

                with open(credentials_path, 'w') as outfile:
                    tomlkit.dump(doc, outfile)

            controller_logger.debug("Streamlit credentials configured successfully")
        except Exception as e:
            controller_logger.error(f"Failed to setup Streamlit credentials: {e}")
            raise
