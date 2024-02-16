import pandas as pd
import numpy as np
import laspy
import os
import subprocess
import re
import json
import stat
from matplotlib.colors import to_hex, to_rgb
from pathlib import Path
import random


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

    out_filename = os.path.join(las_path, 'izzie.las')
    if not os.path.exists(os.path.dirname(out_filename)):
        os.makedirs(os.path.dirname(out_filename))
    las.write(out_filename)
    # print('las file saved at: %s ' % out_filename)


def build_octree(my_path):
    # output_dir = os.path.join(my_path, 'octree', 'Mathieu_z')
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)

    root_dir = os.path.join('..', '..')
    lib = "LD_PRELOAD=%s" % os.path.join(root_dir, 'static', 'PotreeConverter_linux_x64', 'liblaszip.so')
    exe = os.path.join(root_dir, 'static', 'PotreeConverter_linux_x64', 'PotreeConverter')
    input_folder = os.path.join(root_dir, 'src', 'viewer', 'my_test')
    output_folder = os.path.join(root_dir, 'static', 'PotreeConverter_linux_x64')
    opts = " - m  poisson"

    # change file permissions to allow executing
    st = os.stat(exe)
    os.chmod(exe, st.st_mode | stat.S_IEXEC)

    # generate now the octree
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


def gene_settings():
    path_str = os.path.join('..', '..', 'static', '2D', 'viewer', 'js', 'glyphConfig.js')

    # finally get them as a dataframe
    df = parse_js(path_str)

    # add the classification column (ie the glyph id)
    _, classification = np.unique(df.glyphName, return_inverse=True)
    df['classification'] = classification

    # attach now the rgb colors
    out = rgb_helper(df)

    return out


def rgb_helper(df):
    rgb = pd.DataFrame(df.color.map(lambda x: to_rgb(x)).tolist(), columns=['r', 'g', 'b'])
    rgb = 255 * rgb
    rgb = rgb.astype(np.uint8)
    out = pd.concat([df, rgb], axis=1)
    return out


def build_pointcloud(spots_df):
    gs = gene_settings()
    spots_df = spots_df.merge(gs, how='left', left_on='Gene', right_on="gene")

    # fill the nans with the generic values
    generic = gs[gs.gene == 'generic']
    fields = ['color', 'glyphName', 'classification', 'r', 'g', 'b']
    for f in fields:
        spots_df[f] = spots_df[f].fillna(generic[f].values[0])

    spots = spots_df.rename(columns={'Gene_id': 'pointSourceID'})

    target_path = r'my_test'
    build_las(spots, target_path)
    build_octree(target_path)
    print(gs)


def cell_gene_counts(spots_df):
    df = spots_df[['Gene', 'x', 'y', 'z', 'neighbour']]

    df = df[df.neighbour > 0]
    df = df.rename({'Gene': 'gene'}, axis=1)

    output_dir = 'cellData'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for n in np.unique((df.neighbour)):
        temp = df[df.neighbour == n]
        temp.to_json(os.path.join(output_dir,  '%d.json' % n), orient='records')


def cellData_rgb(cellData):
    cellData['Genenames'] = cellData['Genenames'].apply(lambda x: eval(x))
    cellData['CellGeneCount'] = cellData['CellGeneCount'].apply(lambda x: eval(x))
    cellData['ClassName'] = cellData['ClassName'].apply(lambda x: eval(x))
    cellData['Prob'] = cellData['Prob'].apply(lambda x: eval(x))
    cellData['ellipsoid_border'] = cellData['ellipsoid_border'].apply(lambda x: eval(x))
    cellData['sphere_scale'] = cellData['sphere_scale'].apply(lambda x: eval(x))
    cellData['sphere_rotation'] = cellData['sphere_rotation'].apply(lambda x: eval(x))

    cellData['BestClass'] = cellData.apply(lambda r: r['ClassName'][np.argmax(r['Prob'])], axis=1)
    cellData['BestProb'] = cellData.apply(lambda r: np.max(r['Prob']), axis=1)
    cellData = cellData[['Cell_Num', 'X', 'Y', 'Z', 'ClassName', 'Prob', 'sphere_scale', 'sphere_rotation', 'BestClass']]

    uCell_classes = np.unique(sum([d for d in cellData.ClassName], []))
    cellData = cellData.assign(class_id=cellData.ClassName.map(lambda x: [np.where(uCell_classes == d)[0][0] for d in x]))


    path_str = os.path.join('..', '..', 'static', '2D', 'viewer', 'js', 'classConfig.js')
    classConfig = parse_js(path_str)


    # classColours = pd.DataFrame(data)
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

    cellData.to_csv('cellData_rgb.tsv', sep='\t', index=False)



if __name__ == "__main__":
    spots = pd.read_csv(r'/media/dimitris/New Volume/data/Mathieu/WT94_alpha072/pciSeq/data/geneData.tsv', sep='\t')
    build_pointcloud(spots)
    cell_gene_counts(spots)

    # cells = pd.read_csv(r'/media/dimitris/New Volume/data/Mathieu/WT94_alpha072/pciSeq/data/cellData.tsv', sep='\t')
    # cellData_rgb(cells)

    print('Done')