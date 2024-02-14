import pandas as pd
import numpy as np
import laspy
import os
import subprocess
import re
from matplotlib.colors import to_hex, to_rgb
import random


def build_las(data, las_path):
    raw = data.values
    xyz = np.ascontiguousarray(data[['x', 'y', 'z']].values)
    rgb = np.ascontiguousarray(data[['r', 'g', 'b']].values)
    classification = np.ascontiguousarray(data.classification.values, dtype='float32')
    pointSourceID = np.ascontiguousarray(data.pointSourceID.values, dtype='float32')

    # xyz = np.ascontiguousarray(raw[:, 1:4], dtype='float32')
    # rgb = np.ascontiguousarray(raw[:, 7:10], dtype='float32')
    # classification = np.ascontiguousarray(raw[:, 6], dtype='float32')
    # pointSourceID = np.ascontiguousarray(raw[:, -1], dtype='float32')

    hdr = laspy.LasHeader(version="1.4", point_format=7)
    mins = np.floor(np.min(xyz, axis=0))
    np.min(spots[['x', 'y', 'z']].values, axis=0)
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

    # out_path = r"gene_pointclouds_z_spacing_1.5micron/las/%s.las" % gene
    out_filename = os.path.join(las_path, 'izzie.las')
    if not os.path.exists(os.path.dirname(out_filename)):
        os.makedirs(os.path.dirname(out_filename))
    las.write(out_filename)
    # print('las file saved at: %s ' % out_filename)


def build_octree(my_path):
    exe = r"D:\Home\Dimitris\OneDrive - University College London\dev\Javascript\pciSeq_3dviewer\PotreeConverter_windows_x64\PotreeConverter.exe"
    output_dir = os.path.join(my_path, 'octree', 'Mathieu_z')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    result = subprocess.run([exe, my_path, "-o", output_dir, "-m", "poisson"], capture_output=True, shell=True)
    print(result)


def gene_settings():
    path_str = os.path.join('..', '..', 'static', '2D', 'viewer', 'js', 'glyphConfig.js')
    f = open(path_str, "r", encoding='utf8')
    text = f.read()

    # find the text between the square brackets
    sq_brackets = re.findall("\[(.*?)\]", text, flags=re.DOTALL)[0]

    # find now the text between the curly brackets
    crly_brackets = re.findall("{(.*?)}", sq_brackets, flags=re.DOTALL)

    # finally get them as a dataframe
    df = pd.DataFrame([eval("{" + d + "}") for d in crly_brackets])

    # add the pointsource column (ie the gene id)
    # _, pointSourceID = np.unique(df.gene, return_inverse=True)
    # df['pointSourceID'] = pointSourceID

    # add the classification column (ie the glyph id)
    _, classification = np.unique(df.glyphName, return_inverse=True)
    df['classification'] = classification

    # attach now the rgb colors
    rgb = pd.DataFrame(df.color.map(lambda x: to_rgb(x)).tolist(), columns=['r', 'g', 'b'])
    rgb = 255 * rgb
    rgb = rgb.astype(np.uint8)
    out = pd.concat([df, rgb], axis=1)

    return out


def build_pointcloud(spots_df):
    spots = pd.read_csv(r'E:\data\Mathieu\WT94_alpha072\pciSeq\data\geneData.tsv', sep='\t')
    gs = gene_settings()
    spots = spots.merge(gs, how='left', left_on='Gene', right_on="gene")

    # fill the nans with the generic values
    generic = gs[gs.gene == 'generic']
    fields = ['color', 'glyphName', 'classification', 'r', 'g', 'b']
    for f in fields:
        spots[f] = spots[f].fillna(generic[f].values[0])

    spots = spots.rename(columns={'Gene_id': 'pointSourceID'})

    target_path = r'.\test_2'
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


if __name__ == "__main__":
    spots = pd.read_csv(r'E:\data\Mathieu\WT94_alpha072\pciSeq\data\geneData.tsv', sep='\t')
    build_pointcloud(spots)
    cell_gene_counts(spots)