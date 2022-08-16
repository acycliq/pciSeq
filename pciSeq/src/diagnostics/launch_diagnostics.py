import os
import streamlit.bootstrap
from streamlit import config as _config

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'diagnostics.py')

# _config.set_option("server.headless", True)
args = []


streamlit.bootstrap.run(filename, '', args, flag_options={})