import pandas as pd
import numpy as np
import subprocess
from email.parser import BytesHeaderParser
import shutil
import json
import os
import glob
import csv
from pciSeq.src.cell_call.utils import get_out_dir
from pciSeq.src.cell_call.log_config import logger


def make_config_base(dst):
    cellData_tsv = os.path.join(dst, 'data', 'cellData.tsv')
    geneData_tsv = os.path.join(dst, 'data', 'geneData.tsv')

    cellData_dict = {"mediaLink": "../../data/cellData.tsv", "size": str(os.path.getsize(cellData_tsv))}
    geneData_dict = {"mediaLink": "../../data/geneData.tsv", "size": str(os.path.getsize(geneData_tsv))}

    return {
        'cellData': cellData_dict,
        'geneData': geneData_dict,
    }


def make_config_js(dst, w, h):
    appDict = make_config_base(dst)
    cellBoundaries_tsv = os.path.join(dst, 'data', 'cellBoundaries.tsv')
    cellBoundaries_dict = {"mediaLink": "../../data/cellBoundaries.tsv", "size": str(os.path.getsize(cellBoundaries_tsv))}
    roi_dict = {"x0": 0, "x1": w, "y0": 0, "y1": h}
    appDict['cellBoundaries'] = cellBoundaries_dict
    appDict['roi'] = roi_dict
    appDict['zoomLevels'] = 10
    appDict['tiles'] = "https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg"

    config_str = "// NOTES: \n" \
                 "// 1. paths are with respect to the location of 'streaming-tsv-parser.js \n" \
                 "// 2. roi is the image size in pixels. Leave x0 and y0 at zero and set x1 to the width and y1 to the height \n" \
                 "// 3. tiles should point to the folder that keeps your pyramid of tiles. If you do not have that just \n" \
                 "//    change the link to a blind one (change the jpg extension for example). The viewer should work \n" \
                 "//    without the dapi background though \n" \
                 "// 4. size is the tsv size in bytes. I use os.path.getsize() to get it. Not crucial if you \n" \
                 "//    dont get it right, ie the full tsv will still be parsed despite this being wrong. It \n" \
                 "//    is used by the loading page piecharts to calc how far we are \n" \
                 "// 5. Leave zoomLevels to 10 \n" \
                 " function config() { return %s }" % json.dumps(appDict)
    config = os.path.join(dst, 'viewer', 'js', 'config.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    logger.info(' viewer config saved at %s' % config)


def copy_viewer_code(cfg):
    p = subprocess.run(['pip', 'show', 'pciSeq'], stdout=subprocess.PIPE)
    h = BytesHeaderParser().parsebytes(p.stdout)
    pciSeq_dir = os.path.join(h['Location'], 'pciSeq')
    dim = '2D'
    src = os.path.join(pciSeq_dir, 'static', dim)
    dst = get_out_dir(cfg['output_path'], '')

    shutil.copytree(src, dst, dirs_exist_ok=True)
    logger.info(' viewer code (%s) copied from %s to %s' % (dim, src, dst))
    return dst


def splitter_mb(df, dir_path, mb_size):
    """ Splits a text file in (almost) equally sized parts on the disk. Assumes that there is a header in the first line
    :param filepath: The path of the text file to be broken up into smaller files
    :param mb_size: size in MB of each chunk
    :return:
    """
    # OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    else:
        files = glob.glob(dir_path + '/*.*')
        for f in files:
            os.remove(f)

    n = 0
    header_line = df.columns.tolist()
    # header_line = next(handle)[1].tolist()
    file_out, handle_out = _get_file(dir_path, n, header_line)
    # data_row = next(handle)[1].tolist()
    for index, row in df.iterrows():
        row = row.tolist()
        size = os.stat(file_out).st_size
        if size > mb_size*1024*1024:
            logger.info('saved %s with file size %4.3f MB' % (file_out, size/(1024*1024)))
            n += 1
            handle_out.close()
            file_out, handle_out = _get_file(dir_path, n, header_line)
        write = csv.writer(handle_out, delimiter='\t')
        write.writerow(row)

    # print(str(file_out) + " file size = \t" + str(size))
    logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
    handle_out.close()


def splitter_mb(filepath, mb_size):
    """ Splits a text file in (almost) equally sized parts on the disk. Assumes that there is a header in the first line
    :param filepath: The path of the text file to be broken up into smaller files
    :param mb_size: size in MB of each chunk
    :return:
    """
    handle = open(filepath, 'r')
    OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.*')
        for f in files:
            os.remove(f)

    n = 0
    header_line = next(handle)
    file_out, handle_out = _get_file(OUT_DIR, filepath, n, header_line)
    for line in handle:
        size = os.stat(file_out).st_size
        if size > mb_size*1024*1024:
            print('saved %s with file size %4.3f MB' % (file_out, size/(1024*1024)))
            n += 1
            handle_out.close()
            file_out, handle_out = _get_file(OUT_DIR, filepath, n, header_line)
        handle_out.write(str(line))

    # print(str(file_out) + " file size = \t" + str(size))
    print('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
    handle_out.close()


def splitter_n(filepath, n):
    """ Splits a text file in n smaller files
    :param filepath: The path of the text file to be broken up into smaller files
    :param n: determines how many smaller files will be created
    :return:
    """
    filename_ext = os.path.basename(filepath)
    [filename, ext] = filename_ext.split('.')

    OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if ext == 'json':
        df = pd.read_json(filepath)
    elif ext == 'tsv':
        df = pd.read_csv(filepath, sep='\t')
    else:
        df = None

    df_list = np.array_split(df, n)
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.'+ext)
        for f in files:
            os.remove(f)

    for i, d in enumerate(df_list):
        fname = os.path.join(OUT_DIR, filename + '_%d.%s' % (i, ext))
        if ext == 'json':
            d.to_json(fname,  orient='records')
        elif ext == 'tsv':
            d.to_csv(fname, sep='\t', index=False)

