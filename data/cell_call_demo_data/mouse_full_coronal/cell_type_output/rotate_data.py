''''
code used to rotate the data 90 clockwise
'''
import numpy as np
import pandas as pd
import config
import os

cells = pd.read_json('cellData.json')
genes = pd.read_json('geneData.json')

## clockwise where Y is points down.
## From https://math.stackexchange.com/questions/1330161/how-to-rotate-points-through-90-degree
cells_rot = cells.copy()
cells_rot.X = -1 * cells.Y + 27352
cells_rot.Y = cells.X + 7084

## clockwise where Y is points down
## https://math.stackexchange.com/questions/1330161/how-to-rotate-points-through-90-degree
genes_rot = genes.copy()
genes_rot.x = -1 * genes.y + 27352
genes_rot.y = genes.x + 7084


cells_rot.to_json('cellData_rot.json', orient='records')
genes_rot.to_json('geneData_rot.json', orient='records')

coords_path = os.path.join(config.ROOT_DIR, 'dashboard', 'cell_coords.json')
cell_coords = pd.read_json(coords_path)
cell_coords.head()

out = []
for x in cell_coords.itertuples():
    d = np.array(x.coords)
    # print(d)
    if d.size > 1:
        # print(d)
        dr = np.array([-d[:, 1], d[:, 0]]).T
        dr = dr + np.array([27352, 7084])  # because 27352 - 20268 = 7084
        dr = dr.tolist()
        out.append(dr)
    else:
        out.append(None)

cell_coords['coords_rot'] = out

cell_coords_rot = cell_coords.copy()
cell_coords_rot = cell_coords_rot[['cell_id', 'label', 'coords_rot']]
cell_coords_rot = cell_coords_rot.rename(columns={'coords_rot': 'coords'})

cell_coords_rot.to_json('cell_coords_rot.json', orient='records')

print('Done')