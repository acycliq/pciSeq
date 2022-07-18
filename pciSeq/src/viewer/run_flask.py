"""
Function that will open an html page in the browser
"""

from flask import Flask, render_template
import json
import os
import pandas as pd
import numpy as np
import webbrowser
import platform
import random
from threading import Timer
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

dir_path = os.path.dirname(os.path.realpath(__file__))


def get_browser(port_num):
    url = 'http://127.0.0.1:%d' % port_num
    my_os = platform.system()

    if my_os == 'Windows':
        chrome_path = 'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe'
    elif my_os == 'Darwin':  # Is this always the case for MacOS?
        chrome_path = 'open -a /Applications/Google\ Chrome.app'
    elif my_os == 'Linux':
        chrome_path = '/usr/bin/google-chrome'
    else:
        chrome_path = None

    logger.info('Platform is %s' % my_os)
    if chrome_path:
        logger.info('Chrome path: %s' % chrome_path)

    if chrome_path and os.path.isfile(chrome_path):
        webbrowser.register('chrome', None, webbrowser.BackgroundBrowser(chrome_path),  preferred=True)
        wb = webbrowser.get('chrome').open_new_tab(url)
    else:
        wb = webbrowser.open_new(url)

    if not wb:
        logger.info('Could not open browser')


def open_browser():
    webbrowser.open_new('http://127.0.0.1:5000/')


def flask_app_start(config):
    template_dir = os.path.abspath('../../../viewer')
    port = 5000 + random.randint(0, 999)
    # viewer / genes_datatable.html
    # pciSeq/src/viewer/run_flask.py
    flask_app = Flask(__name__,
                static_url_path='',         # remove the static folder path
                static_folder=r'..\..\..\viewer',   # set here the path of the folder to be served. The js files referenced in your html file are with respect to this folder. Adjust the paths in your html file (look for the <script src="some/path/file.js"></script> tag, so that the js libraries will be parsed
                template_folder=r'D:\Home\Dimitris\OneDrive - University College London\dev\Python\pciSeq\viewer')  # set here the path to the folder where your html page lives. Absolute and relative paths both work fine

    @flask_app.route("/")
    def index():
        return render_template("genes_datatable.html", data=None)

    Timer(1, get_browser, [port]).start()
    flask_app.run(port=port)


def mk_ini(cellData, geneData):
    ini = {}
    # ini['name'] = json.dumps(config.DEFAULT['dataset_id'])
    # ini['roi'] = json.dumps(config.DEFAULT['roi'])
    # ini['imageSize'] = json.dumps(pixel_dims)
    ini['cellData'] = cellData.to_json(orient='records')
    ini['geneData'] = geneData.to_json(orient='records')
    # ini['tiles'] = json.dumps(os.path.join(tiles_root_path, '{z}', '{y}', '{x}.png'))

    return ini

if __name__ == "__main__":
    geneData = pd.read_csv(r"C:\Users\skgtdni\AppData\Local\Temp\pciSeq\geneData.tsv", sep='\t')
    cellData = pd.read_csv(r"C:\Users\skgtdni\AppData\Local\Temp\pciSeq\cellData.tsv", sep='\t')
    ini = mk_ini(cellData, geneData)
    flask_app_start(ini)
