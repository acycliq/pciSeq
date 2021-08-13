"""
script that reads from dropbox (using dropbox API) cell_by_gene.csv, manifest.json, micron_to_mosaic_pixel_transform.csv
and saves them locally on the disk
"""
import dropbox
import pathlib
import pandas as pd
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def get_dbx_paths(dbx, ROOT_PATH, cursor=False):
    res = []
    if cursor is False:
        out = dbx.files_list_folder(ROOT_PATH, recursive=True)
    else:
        out = dbx.files_list_folder_continue(cursor)
    for p in out.entries:
        filename = p.path_display
        ext = pathlib.Path(filename).suffix
        visited = pd.read_csv('dbx_fetch.log').values
        if filename not in visited:
            if (ext in ['.csv', '.json']) and (os.path.basename(filename) != 'detected_transcripts.csv') and (filename not in visited):
                download_file(dbx, filename)
                with open('dbx_fetch.log', 'a') as fd:
                    fd.write(filename)
                    fd.write("\n")
                res.append(filename)
        else:
            logger.info('Skipping. File: %s already downloaded' % filename)
    if out.has_more:
        get_dbx_paths(dbx, ROOT_PATH, cursor=out.cursor)
    return res


def download_file(dbx, filename):
    fName = filename.replace('/MERFISH:DRONC', 'D:/vizgen_raw_data')
    parent_dir = Path(fName).parent.absolute()
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)
    with open(fName, "wb") as f:
        logger.info('Downloading %s' % filename)
        metadata, res = dbx.files_download(path=filename)
        f.write(res.content)
        logger.info('Saved %s' % fName)


def worker(merfish_id, dbx):
    pid = os.getpid()
    print(f'Worker pid: {pid}, processing {merfish_id} ')
    DBX_ROOT_PATH = '/MERFISH:DRONC/' + merfish_id
    get_dbx_paths(dbx, DBX_ROOT_PATH)


def run_par(dbx, merfish_ids):
    """
    same as run() but uses multiprocessing
    """
    processes = cpu_count()
    with Pool(processes=processes) as pool:
        pool.map(partial(worker, dbx=dbx), merfish_ids)


def run(dbx, merfish_ids):
    for merfish_id in merfish_ids:
        DBX_ROOT_PATH = '/MERFISH:DRONC/' + merfish_id
        get_dbx_paths(dbx, DBX_ROOT_PATH)


if __name__ == "__main__":
    DROPBOX_TOKEN = "Xlk6baOMGrMAAAAAAAAAATuJFoFOL2_2ix6sa3dDK2IuQRarxEQD3EG9SMIig1FS"
    merfish_ids = ['MERFISH_F_E', 'MERFISH_F_F', 'MERFISH_M_C', 'MERFISH_M_Z']
    dbx = dropbox.Dropbox(DROPBOX_TOKEN)

    run(dbx, merfish_ids)
    print('Done!')

