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
cells_landscape = cells.copy()
cells_landscape.X = cells.Y - 27352/2 + 20268/2 + 7084/2
cells_landscape.Y = -(cells.X - 20268/2) + 27352/2 + 7084/2

## clockwise where Y is points down
## https://math.stackexchange.com/questions/1330161/how-to-rotate-points-through-90-degree
genes_landscape = genes.copy()
genes_landscape.x = genes.y - 27352/2 + 20268/2 + 7084/2
genes_landscape.y = -(genes.x - 20268/2) + 27352/2 + 7084/2


cells_landscape.to_json('cellData_landscape.json', orient='records')
genes_landscape.to_json('geneData_landscape.json', orient='records')

coords_path = os.path.join(config.ROOT_DIR, 'dashboard', 'cell_coords.json')
cell_coords = pd.read_json(coords_path)
cell_coords.head()

out = []
for x in cell_coords.itertuples():
    d = np.array(x.coords)
    # print(d)
    if d.size > 1:
        # print(d)
        dr = np.array([d[:, 1], -d[:, 0]]).T
        dr = dr + np.array([(20268 - 27352) / 2, (20268 + 27352) / 2])
        dr = dr + np.array([7084/2, 7084/2])  # because 27352 - 20268 = 7084
        dr = dr.astype(np.int).tolist()
        out.append(dr)
    else:
        out.append(None)

cell_coords['coords_landscape'] = out

cell_coords_landscape = cell_coords.copy()
cell_coords_landscape = cell_coords_landscape[['cell_id', 'label', 'coords_landscape']]
cell_coords_landscape = cell_coords_landscape.rename(columns={'coords_landscape': 'coords'})
cell_coords_landscape.cell_id = cell_coords_landscape.cell_id.astype(np.int)
cell_coords_landscape.label = cell_coords_landscape.label.astype(np.int)

cell_coords_landscape.to_json('cell_coords_landscape.json', orient='records')

print('Done')