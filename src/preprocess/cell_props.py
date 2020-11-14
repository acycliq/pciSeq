import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import skimage.measure as skmeas
from collections import defaultdict
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()


def calc_props(stage):
    """ For the clipped cells, collate together all the necessary label_image arrays so that the cell is
    not clipped anymore. Then read that mosaic of arrays and get cell centroid, area and other properties.
    For the unclipped cells these properties are kept in image_objects.
    Combine all (clipped and unclipped) and make a single dataframe with columns: 'label', 'tile_id', 'area',
    'x_local', 'y_local', 'tile_offset_x', 'tile_offset_y', 'is_clipped', 'x', 'y', 'coords'.
    It also assigns a unique global label to the cells. The last column with label 'coords' keeps the
    coordinates of the cell boundary
    """

    # make a dict where each key is a list of tiles and each value is a list of labels. The latter
    # holds the labels whose corresponding objects as these were identified by the cell segmentation
    # extend outside one single tile and span across the tiles kept in the corresponding key of the dict
    v = defaultdict(list)
    for key, value in sorted(list(stage.merge_register.entries.items())):
        if key not in v[tuple(value)]:
            v[tuple(value)].append(key)

    # first, do the merge cells
    # check the register to see which cells are merged
    merged_props_dict = {'label': [], 'tile_id': [], 'area': [], 'x_local': [], 'y_local': []}
    for tt in v.items():
        labels = tt[1]
        tile_ids = sorted(tt[0])

        logger.info('calculating centroid and area for the merged cells with labels %s' % labels)
        if len(tile_ids) > 1:
            label_image = stage.collate_arrays(tile_ids)
        else:
            label_image = None  # you shouldnt get here. Since we merge cells we need at least two tiles

        # sanity check
        coo_mat = set(coo_matrix(label_image).data)
        assert np.all([label in coo_mat for label in labels]), 'a label that should exist in the collated label_image seem to be missing'

        props = skmeas.regionprops(label_image.astype(np.int32))
        clipped_cells = list(filter(lambda d: d.label in labels, props))
        # unclipped_cells = list(filter(lambda d: d.label not in labels, props))
        for i, p in enumerate(clipped_cells):
            centroid_tile_id, centroid_coords = locate_tile(stage, p.centroid, tile_ids)

            logger.info('cell with label %d is clipped by tile_ids %s' % (p.label, tile_ids))

            # # sanity check
            # filtered_tiles = list(filter(lambda d: d['tile_id'] in tile_ids, self.tiles))
            # mask = [d['image_objects'].label == labels[i] for d in filtered_tiles]
            # total_area = sum([d['image_objects'].area[mask[i]].values[0] for i, d in enumerate(filtered_tiles)])
            # # logger.info('total sum for the sanity check is: %d' % total_area)
            # assert int(total_area) == int(p.area)

            merged_props_dict['label'].append(p.label)
            merged_props_dict['tile_id'].append(centroid_tile_id)
            merged_props_dict['area'].append(p.area)
            merged_props_dict['x_local'].append(centroid_coords[1])
            merged_props_dict['y_local'].append(centroid_coords[0])

            # sort the dict by the label
            merged_props_dict = pd.DataFrame(merged_props_dict).sort_values(by=['label']).to_dict(orient='list')
            # logger.info('merged_props_dict ends')

    img_obj = pd.concat([d['image_objects'] for d in stage.tiles])
    dup_vals = list(set(img_obj[img_obj.duplicated('label')].label.values))
    # np.setdiff1d(np.array(self.cell_props['label']), dup_vals) # <--  WHY I HAVE THAT?

    unclipped_labels = list(set(img_obj.label.values) - set(merged_props_dict['label']))
    unclipped_cells = img_obj.iloc[np.isin(img_obj.label.values, unclipped_labels), :]

    # sanity check. Cannot have duplicate labels. If a label is duplicate then the cell is clipped, hence it should have been captured by
    # 'merged_props_dict'
    assert unclipped_cells[unclipped_cells.duplicated('label')].empty, "Dataframe 'unclipped_img_obj' contains dublicate labels."

    # unclipped_cells = img_obj.drop_duplicates('label', keep=False)
    merged = pd.DataFrame(merged_props_dict)

    df_1 = unclipped_cells.merge(pd.DataFrame(stage.tiles)[['tile_id', 'tile_offset_x', 'tile_offset_y']], how='left', on='tile_id')
    df_1['is_clipped'] = False
    df_2 = merged.merge(pd.DataFrame(stage.tiles)[['tile_id', 'tile_offset_x', 'tile_offset_y']], how='left', on='tile_id')
    df_2['is_clipped'] = True

    cell_props = pd.concat([df_1.reset_index(drop=True), df_2.reset_index(drop=True)], axis=0).reset_index(drop=True)
    cell_props['x'] = cell_props.x_local + cell_props.tile_offset_x
    cell_props['y'] = cell_props.y_local + cell_props.tile_offset_y

    # sort by x and y and then relabel so that labelling looks more deterministic.
    # Not necessary, it is nice to have tho. If however that takes place then the
    # label_image array has to be updated to be inline with the new labels
    cell_props = cell_props.sort_values(by=['x', 'y'])

    global_label = np.arange(cell_props.shape[0]) + 1
    label_map = np.column_stack((cell_props.label.values, global_label))
    cell_props['label'] = global_label
    cell_props['label'] = global_label

    logger.info('')
    logger.info('updating the labels on label_image to align them with cell_props')
    for tile in stage.tiles:
        data = update_label_image(tile, label_map)
        tile['label_image'].data = data.astype(np.int64)

    logger.info('Done')

    # relabel now the key in the merge_register dict to align them
    # with the global labels
    stage.merge_register = _update_dict_keys(stage, label_map)

    # DO I ALSO NEED TO UPDATE THE IMAGE_OBJECTS FOR EACH FOV??
    # Labels there need to be updated for consistency
    # However in a cell is split in to two parts and exists in two
    # tiles then x_local and y_local coords show the centroid of the partial shape
    # The array cell_props is the correct place for such lookups.

    # Find now the outline of the cells
    _cell_boundaries = stage.cell_boundaries(cell_props)
    out = cell_props.merge(_cell_boundaries, how='left', on=['label'])

    return out


def locate_tile(stage, centroid, tile_ids):
    # first get the tile_id
    t = stage.tile_topo(tile_ids)
    row = int(centroid[0] // stage.tile_shape[0])  # <-- I think I should be dividing by tile_shape[1] instead
    col = int(centroid[1] // stage.tile_shape[1])  # <-- I think I should be dividing by tile_shape[0] instead

    tile_id = t[row, col]
    assert ~np.isnan(tile_id)
    assert tile_id in tile_ids

    # calc the local coordinates
    coord_row = centroid[0] % stage.tile_shape[0]  # <-- I think I should be dividing by tile_shape[1] instead
    coord_col = centroid[1] % stage.tile_shape[1]  # <-- I think I should be dividing by tile_shape[0] instead

    return tile_id, (coord_row, coord_col)


def update_label_image(tile, label_map):
    _x = tile['label_image'].data
    _y = label_map[:, 0]
    idx = _remap(_x, _y)
    return label_map[idx, 1]


def _remap(x, y):
    index = np.argsort(y)
    sorted_y = y[index]
    sorted_index = np.searchsorted(sorted_y, x)

    xindex = np.take(index, sorted_index, mode="clip")
    mask = y[xindex] != x

    result = np.ma.array(xindex, mask=mask)
    return result


def _update_dict_keys(stage, label_map):
    # relabel now the key in the merge_register dict to align them
    # with the global labels
    dict_keys = np.array(list(stage.merge_register.entries.keys()))
    keys_mask = _remap(dict_keys, label_map[:, 0])
    global_keys = label_map[keys_mask, 1]
    dict_vals = list(stage.merge_register.entries.values())
    global_dict = dict(zip(global_keys, dict_vals))
    return global_dict

