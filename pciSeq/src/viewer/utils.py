import pandas as pd
import numpy as np
from pciSeq.src.core.utils import get_pciSeq_install_dir
import shutil
import json
import os
import glob
import csv
from pathlib import Path
from matplotlib.colors import to_hex, to_rgb
import re
import laspy
import subprocess
import stat
import logging

viewer_utils_logger = logging.getLogger(__name__)


def make_config_base(dst):
    cellData_tsv = os.path.join(dst, 'data', 'cellData.tsv')
    geneData_tsv = os.path.join(dst, 'data', 'geneData.tsv')

    cellData_dict = {"mediaLink": "../../data/cellData.tsv", "size": str(os.path.getsize(cellData_tsv))}
    geneData_dict = {"mediaLink": "../../data/geneData.tsv", "size": str(os.path.getsize(geneData_tsv))}

    return {
        'cellData': cellData_dict,
        'geneData': geneData_dict,
    }


def make_config_js(dst, img_shape):
    appDict = make_config_base(dst)
    cellBoundaries_tsv = os.path.join(dst, 'data', 'cellBoundaries.tsv')
    cellBoundaries_dict = {"mediaLink": "../../data/cellBoundaries.tsv",
                           "size": str(os.path.getsize(cellBoundaries_tsv))}
    roi_dict = {"x0": 0, "x1": img_shape[1], "y0": 0, "y1": img_shape[0]}
    appDict['cellBoundaries'] = cellBoundaries_dict
    appDict['roi'] = roi_dict
    appDict['maxZoom'] = 8
    appDict['layers'] = {
        'empty': "",
        'dapi': "https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg",

    }
    appDict['spotSize'] = 1 / 16

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
    config = os.path.join(dst, 'viewer', 'libs', 'config.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    viewer_utils_logger.info('viewer config saved at %s' % config)


def make_classConfig_nsc_js(labels, dst):
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
    config = os.path.join(dst, 'viewer', 'libs', 'classConfig.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    viewer_utils_logger.info('cell class color scheme saved at %s' % config)


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
    viewer_utils_logger.info(' classConfig saved at %s' % config)


def make_classConfig_js(pciSeq_dir, dst):
    json_file = os.path.join(pciSeq_dir, 'static', 'color_scheme', 'classColors.json')
    with open(json_file, 'r') as f:
        data = json.load(f)

    config_dict = [d for d in data if '//' not in d.keys()]
    config_str = " function classColorsCodes() { return %s }" % json.dumps(config_dict)
    config = os.path.join(dst, 'viewer', 'libs', 'classConfig.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    viewer_utils_logger.info('cell class color scheme saved at %s' % config)


def make_glyphConfig_js(gene_panel, pciSeq_dir, dst):
    json_file = os.path.join(pciSeq_dir, 'static', 'color_scheme', 'geneColors.json')
    with open(json_file, 'r') as f:
        data = json.load(f)
        out = _make_glyphConfig_js(gene_panel, data, dst)
    return out


def _make_glyphConfig_js(gene_panel, data, dst):

    config_dict = [d for d in data if '//' not in d.keys()]

    # read the preset settings (gene color and glyph)
    df = pd.DataFrame(config_dict)

    # for your gene panel get the color/glyph for the preset file
    in_gene_panel = np.in1d(df.gene.values, gene_panel)
    df1 = df[in_gene_panel]

    # use the generic settings for those genes in the gene panel that do not have color/glyph,
    unspecified = set(gene_panel) - set(df.gene)
    generic = df[df.gene == 'generic']
    df2 = pd.DataFrame({'gene': list(unspecified),
                  'color': generic.color.tolist() * len(unspecified),
                  'glyphName': generic.glyphName.tolist() * len(unspecified),
                  })

    # keep now everything in one dataframe
    config_df = pd.concat([df1, df2], ignore_index=True)

    # and save the color scheme to the glyphConfig.js file
    config_str = " function glyphSettings() { return %s }" % json.dumps(config_df.to_dict(orient='records'))
    config = os.path.join(dst, 'viewer', 'libs', 'glyphConfig.js')
    with open(config, 'w') as data:
        data.write(str(config_str))
    viewer_utils_logger.info('glyph color scheme saved at %s' % config)
    return data


def copy_viewer_code(dst, dim):
    pciSeq_dir = get_pciSeq_install_dir()
    src = os.path.join(pciSeq_dir, 'static', dim)

    shutil.copytree(src, dst, dirs_exist_ok=True)
    viewer_utils_logger.info('viewer code (%s) copied from %s to %s' % (dim, src, dst))

    return pciSeq_dir


def build_pointcloud(spots_df, pciSeq_dir, dst, gs=None):
    """
    USE FUNCTION OVERLOADING!!! (multidispatch maybe or typing)
    """
    if gs is None:
        gs = gene_settings(pciSeq_dir)
    else:
        gs = _gene_settings(gs)
    spots_df = spots_df.merge(gs, how='left', left_on='Gene', right_on="gene")
    spots_df = spots_df.dropna()
    # fill the nans with the generic values
    generic = gs[gs.gene == 'generic']
    fields = ['color', 'glyphName', 'classification', 'r', 'g', 'b']
    for f in fields:
        spots_df[f] = spots_df[f].fillna(generic[f].values[0])

    spots = spots_df.rename(columns={'Gene_id': 'pointSourceID'})

    data_folder = os.path.join(dst, 'data')
    Path(data_folder).mkdir(parents=True, exist_ok=True)

    build_las(spots, data_folder)
    build_octree(pciSeq_dir, data_folder)
    viewer_utils_logger.info('octree saved at: %s' % data_folder)


def gene_settings(pciSeq_dir):
    json_file = os.path.join(pciSeq_dir, 'static', 'color_scheme', 'geneColors.json')
    with open(json_file, 'r') as f:
        data = json.load(f)

    out = _gene_settings(data)
    return out


def _gene_settings(data):
    config_dict = [d for d in data if '//' not in d.keys()]
    df = pd.DataFrame(config_dict)

    # path_str = os.path.join('pciSeq', 'static', '2D', 'viewer', 'libs', 'glyphConfig.js')

    # finally get them as a dataframe
    # df_2 = parse_js(path_str)

    # add the classification column (ie the glyph id)
    _, classification = np.unique(df.glyphName, return_inverse=True)
    df['classification'] = classification

    # attach now the rgb colors
    return rgb_helper(df)

def build_las(data, las_path):
    raw = data.values
    xyz = np.ascontiguousarray(data[['x', 'y', 'z']].values)
    rgb = np.ascontiguousarray(data[['r', 'g', 'b']].values)
    classification = np.ascontiguousarray(data.classification.values, dtype='float32')
    pointSourceID = np.ascontiguousarray(data.pointSourceID.values, dtype='float32')

    hdr = laspy.LasHeader(version="1.4", point_format=7)
    mins = np.floor(np.min(xyz, axis=0))
    np.min(data[['x', 'y', 'z']].values, axis=0)
    # mins = [352, 6126, 0]
    hdr.offset = mins
    hdr.scales = np.array([0.001, 0.001, 0.001])

    # 2. Create a Las
    las = laspy.LasData(hdr)

    las.x = xyz[:, 0]
    las.y = xyz[:, 1]
    las.z = xyz[:, 2]
    las.red = rgb[:, 0]
    las.green = rgb[:, 1]
    las.blue = rgb[:, 2]
    las.classification = classification
    las.pt_src_id = pointSourceID
    # las.intensity = i

    out_filename = os.path.join(las_path, 'pciSeq.las')
    if not os.path.exists(os.path.dirname(out_filename)):
        os.makedirs(os.path.dirname(out_filename))
    las.write(out_filename)
    # print('las file saved at: %s ' % out_filename)


def build_octree(pciSeq_dir, input_folder):
    # output_dir = os.path.join(my_path, 'octree', 'Mathieu_z')
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)

    lib = "LD_PRELOAD=%s" % os.path.join(pciSeq_dir, 'static', 'PotreeConverter_linux_x64', 'liblaszip.so')
    exe = os.path.join(pciSeq_dir, 'static', 'PotreeConverter_linux_x64', 'PotreeConverter')
    output_folder = os.path.join(input_folder, 'pointclouds')
    opts = " - m  poisson"

    # change file permissions to allow executing
    st = os.stat(exe)
    os.chmod(exe, st.st_mode | stat.S_IEXEC)

    # generate now the octree
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    cmd = [lib + " " + exe + " " + input_folder + " -o " + output_folder + opts]
    result = subprocess.run(cmd, capture_output=True, shell=True)
    print(result)


def parse_js(path_str):
    f = open(path_str, "r", encoding='utf8')
    text = f.read()

    # find the text between the square brackets
    sq_brackets = re.findall("\[(.*?)\]", text, flags=re.DOTALL)[0]

    # find now the text between the curly brackets
    crly_brackets = re.findall("{(.*?)}", sq_brackets, flags=re.DOTALL)

    # finally get them as a dataframe
    try:
        df = pd.DataFrame([eval("{" + d + "}") for d in crly_brackets])
    except:
        df = pd.DataFrame(["{" + d + "}" for d in crly_brackets])
    return df


def rgb_helper(df):
    rgb = pd.DataFrame(df.color.map(lambda x: to_rgb(x)).tolist(), columns=['r', 'g', 'b'])
    rgb = 255 * rgb
    rgb = rgb.astype(np.uint8)
    out = pd.concat([df, rgb], axis=1)
    return out


def cellData_rgb(cellData, pciSeq_dir, dst):
    '''
    Better move this on the javascript side!
    '''
    cellData['BestClass'] = cellData.apply(lambda r: r['ClassName'][np.argmax(r['Prob'])], axis=1)
    cellData['BestProb'] = cellData.apply(lambda r: np.max(r['Prob']), axis=1)

    uCell_classes = np.unique(sum([d for d in cellData.ClassName], []))
    cellData = cellData.assign(class_id=cellData.ClassName.map(lambda x: [np.where(uCell_classes == d)[0][0] for d in x]))

    # path_str = os.path.join(pciSeq_dir, 'static', '2D', 'viewer', 'libs', 'classConfig.js')
    json_file = os.path.join(pciSeq_dir, 'static', 'color_scheme', 'classColors.json')
    with open(json_file, 'r') as f:
        data = json.load(f)

    config_dict = [d for d in data if '//' not in d.keys()]
    classConfig = pd.DataFrame(config_dict)

    # classConfig = parse_js(path_str)

    classConfig = rgb_helper(classConfig)
    cellData = cellData.merge(classConfig, how='left', left_on='BestClass', right_on="className")
    cellData = cellData[['Cell_Num', 'X', 'Y', 'Z', 'ClassName', 'class_id', 'Prob', 'sphere_scale', 'sphere_rotation', 'r', 'g', 'b']]

    cellData = cellData.rename(columns={"Cell_Num": "label",
                                        "X": "x",
                                        "Y": "y",
                                        "Z": "z",
                                        "Prob": "class_prob"})

    # fill the nans with the generic values
    generic = classConfig[classConfig.className == 'Generic']
    fields = ['r', 'g', 'b']
    for f in fields:
        cellData[f] = cellData[f].fillna(generic[f].values[0])

    data_folder = os.path.join(dst, 'data')
    Path(data_folder).mkdir(parents=True, exist_ok=True)
    cellData.to_csv(os.path.join(data_folder, 'cellData_rgb.tsv'), sep='\t', index=False)


def cell_gene_counts(spots_df, dst):
    df = spots_df[['Gene', 'x', 'y', 'z', 'neighbour']]

    df = df[df.neighbour > 0]
    df = df.rename({'Gene': 'gene'}, axis=1)

    data_folder = os.path.join(dst, 'data', 'cell_gene_counts')
    Path(data_folder).mkdir(parents=True, exist_ok=True)
    for n in np.unique((df.neighbour)):
        temp = df[df.neighbour == n]
        temp.to_json(os.path.join(data_folder,  '%d.json' % n), orient='records')

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
            viewer_utils_logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
            n += 1
            handle_out.close()
            file_out, handle_out = _get_file(dir_path, n, header_line)
        write = csv.writer(handle_out, delimiter='\t')
        write.writerow(row)

    # print(str(file_out) + " file size = \t" + str(size))
    viewer_utils_logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
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
