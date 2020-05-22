''''
produces the fov boundaries in geojson format. Writes the fov_boundaries.js which is imported by the index.html
I think I will write something to do that inside jabascript anyway....
'''

import logging
import config
import os
import json
from src.preprocess.fov import Fov

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

dir_path = os.path.dirname(os.path.realpath(__file__))


def side_length(coords):
    x_set = set([d[0] for d in coords])
    y_set = set([d[1] for d in coords])

    assert len(x_set) == len(y_set) == 2

    x_list = sorted(list(x_set))
    y_list = sorted(list(y_set))

    x_side = [d - x_list[i - 1] for i, d in enumerate(x_list) if i > 0]
    y_side = [d - y_list[i - 1] for i, d in enumerate(y_list) if i > 0]

    assert x_side == y_side

    return int(x_side[0])


def fov_bounding_box(fovs):
    '''
    makes the polygons for the blocks.
    :param my_config:
    :return:
    '''

    bbox = [list(zip(d['fov_range']['x'], d['fov_range']['y'])) for d in fovs.fovs]

    coords = [[
                [d[0][0], d[0][1]],
                [d[0][0], d[1][1]],
                [d[1][0], d[1][1]],
                [d[1][0], d[0][1]],
                [d[0][0], d[0][1]]
               ]
               for d in bbox]

    fov_name = [d['fov_id'] for d in fovs.fovs]
    return fov_name, coords


def fov_boundaries(fovs):
    '''
    writes a geoJson file describing the polygons that define each one of the blocks
    :param my_config:
    :return:
    '''

    fov_boundaries_path = os.path.join(config.ROOT_DIR, 'dashboard', 'js', 'fov_boundaries.js')
    fov_name, coords = fov_bounding_box(fovs)

    geometry = []
    for i in range(len(coords)):
        item = {"type": "Polygon"}
        item["coordinates"] = [coords[i]]
        geometry.append(item)

    properties = []
    for i in range(len(coords)):
        item = {'id': fov_name[i]}
        item['side_length'] = side_length(coords[i])
        item['feature_num'] = i
        properties.append(item)

    features = []
    for i in range(len(geometry)):
        item = {"type": "Feature"}
        item["properties"] = properties[i]
        item["geometry"] = geometry[i]
        features.append(item)

    featureCollection = {}
    featureCollection["type"] = "FeatureCollection"
    featureCollection["features"] = features

    jsonData = json.dumps(featureCollection)
    # output = cfg['BLOCK_BOUNDARIES']
    with open(fov_boundaries_path, "w") as text_file:
        print(" var fov_boundaries = %s" % jsonData, file=text_file)

    logger.info('Fov bounding box boundaries saved at %s' % fov_boundaries_path)
    return True


if __name__ == "__main__":
    cfg = config.PREPROCESSOR

    fovs_across = cfg['FOVS_ACROSS']
    fovs_down = cfg['FOVS_DOWN']
    fovs = Fov(fovs_across, fovs_down, cfg)

    fov_boundaries(fovs)
