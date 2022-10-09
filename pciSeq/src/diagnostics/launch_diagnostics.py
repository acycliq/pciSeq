import os
import sys
import subprocess


def launch_dashboard():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')
    exe = sys.executable

    subprocess.Popen([exe, "-m" "streamlit", "run", filename])


if __name__ == "__main__":
    launch_dashboard()
