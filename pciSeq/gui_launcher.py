#!/usr/bin/env python
"""
pciSeq GUI Launcher

This is the main entry point for non-technical users.
After installing pciSeq, simply type 'pciSeq' in your terminal to launch the GUI.

The GUI will open in your browser where you can:
1. Select your data files
2. Configure parameters
3. Start the analysis

The server will continue running until you close this window (Ctrl+C).
"""

import sys
import signal
from pciSeq.src.realtime_viewer.server import RealtimeViewerServer
from pciSeq.src.core.logger import setup_logger


def signal_handler(sig, frame):
    """Handle Ctrl+C gracefully."""
    print("\n\nShutting down pciSeq GUI...")
    sys.exit(0)


def main():
    """Launch the pciSeq GUI server."""
    # Set up logger
    setup_logger()

    print("=" * 60)
    print("pciSeq GUI Launcher")
    print("=" * 60)
    print()
    print("Starting server...")
    print("The GUI will open automatically in your browser.")
    print()
    print("To stop the server, press Ctrl+C")
    print("=" * 60)
    print()

    # Register Ctrl+C handler
    signal.signal(signal.SIGINT, signal_handler)

    # Create and start the server
    # auto_open_browser=True will automatically open the browser
    server = RealtimeViewerServer(port=5001, host="127.0.0.1", auto_open_browser=True)

    server.start()

    print("\nServer started successfully!")
    print("If the browser didn't open automatically, visit: http://127.0.0.1:5001")
    print()

    # Keep the server running
    try:
        # Wait indefinitely (server runs in background thread)
        signal.pause()
    except AttributeError:
        # signal.pause() not available on Windows
        import time

        while True:
            time.sleep(1)


if __name__ == "__main__":
    main()
