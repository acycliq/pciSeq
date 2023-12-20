import pandas as pd
import numpy as np
from pciSeq.src.core.utils import get_pciSeq_install_dir
import shutil
import json
import os
import glob
import csv
from pciSeq.src.core.log_config import logger


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
    cellBoundaries_dict = {"mediaLink": "../../data/cellBoundaries.tsv",
                           "size": str(os.path.getsize(cellBoundaries_tsv))}
    roi_dict = {"x0": 0, "x1": w, "y0": 0, "y1": h}
    appDict['cellBoundaries'] = cellBoundaries_dict
    appDict['roi'] = roi_dict
    appDict['maxZoom'] = 8
    appDict['layers'] = {
        'empty': "",
        'dapi': "https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg",

    }
    appDict['spotSize'] = 1/16

    config_str = "// NOTES: \n" \
                 "// 1. paths in 'cellData', 'geneData' and 'cellBoundaries' are with respect to the location of \n" \
                 "//    'streaming-tsv-parser.js' \n" \
                 "// 2. size is the tsv size in bytes. I use os.path.getsize() to get it. Not crucial if you \n" \
                 "//    don't get it right, ie the full tsv will still be parsed despite this being wrong. It \n" \
                 "//    is used by the loading page piecharts to calc how far we are. \n" \
                 "// 3. roi is the image size in pixels. Leave x0 and y0 at zero and set x1 to the width and y1 to the height. \n" \
                 "// 4. layers is a dict. Each key/value pair contains the string (the name) of the background image and the \n" \
                 "//    location of the folder that the corresponding pyramid of tiles. If the tiles are stored locally, they \n" \
                 "//    should be kept in a folder which is served, for example next to the tsv flatfiles. The path should be \n" \
                 "//    in relation to the location of the index.html If you do not have a pyramid of tiles just \n" \
                 "//    change the link to a blind one (change the jpg extension for example or just use an empty string). \n" \
                 "//    The viewer should work without the dapi background though. \n" \
                 "//    If the dict has more than one entries then a small control with radio button will appear at the top \n" \
                 "//    right of the viewer to switch between different background images. \n" \
                 "// 5. maxZoom: maximum zoom levels. In most cases a value of 8 if good enough. If you have a big image, like \n" \
                 "//    full coronal section for example then a value of 10 would make sense. Note that this should be typically \n" \
                 "//    inline with the zoom level you used when you did \n" \
                 "//    the pyramid of tiles. No harm is it less. If it is greater, then for these extra zoom levels there will \n" \
                 "//    be no background image. \n" \
                 "// 6. spotSize: Scalar. Use this to adjust the screen-size of your spots before they morph into glyphs. \n" \
                 " function config() { return %s }" % json.dumps(appDict)
    config = os.path.join(dst, 'viewer', 'js', 'config.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    logger.info(' viewer config saved at %s' % config)


def make_classConfig_js(labels, dst):
    # remove Zero. It is appended later on
    if 'Zero' in labels:
        labels.remove('Zero')

    colours = ["#f3c300", "#875692", "#f38400", "#a1caf1", "#be0032",
               "#c2b280", "#848482", "#008856", "#e68fac", "#0067a5",
               "#f99379", "#604e97", "#f6a600", "#b3446c", "#dcd300",
               "#882d17", "#8db600", "#654522", "#e25822", "#2b3d26"]
    n = len(colours)
    config_dict = [{'className': labels[i],
                   'IdentifiedType': labels[i],
                   'color': colours[i % n]}
                  for i, v in enumerate(labels)]
    config_dict.append({'className': 'Zero', 'IdentifiedType': 'Zero', 'color': '#000000'})
    config_dict.append({'className': 'Other', 'IdentifiedType': 'Other', 'color': '#C0C0C0'})
    config_str = " function classColorsCodes() { return %s }" % json.dumps(config_dict)
    config = os.path.join(dst, 'viewer', 'js', 'classConfig.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    logger.info(' classConfig saved at %s' % config)


def copy_viewer_code(cfg, dst):
    pciSeq_dir = get_pciSeq_install_dir()
    dim = '2D'
    src = os.path.join(pciSeq_dir, 'static', dim)

    shutil.copytree(src, dst, dirs_exist_ok=True)
    logger.info(' viewer code (%s) copied from %s to %s' % (dim, src, dst))
    return dst


def _get_file(OUT_DIR, filepath, n, header_line):
    [filename, ext] = os.path.basename(filepath).split('.')
    file = os.path.join(OUT_DIR, filename + '_%d.%s' % (n, ext))
    handle = open(file, "a")
    handle.write(header_line)
    return file, handle


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
        if size > mb_size * 1024 * 1024:
            logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
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
        if size > mb_size * 1024 * 1024:
            print('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
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
        files = glob.glob(OUT_DIR + '/*.' + ext)
        for f in files:
            os.remove(f)

    for i, d in enumerate(df_list):
        fname = os.path.join(OUT_DIR, filename + '_%d.%s' % (i, ext))
        if ext == 'json':
            d.to_json(fname, orient='records')
        elif ext == 'tsv':
            d.to_csv(fname, sep='\t', index=False)
