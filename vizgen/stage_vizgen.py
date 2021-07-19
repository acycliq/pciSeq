"""
stages the vizgen data relevant to the cells. It actually reads the rotated cellprops and cell boundaties and keeps
only what is inside the ROI. Then it labels the cells. Label 0 is for the background
"""
import pandas as pd
import numpy as np
import numba



@numba.jit(nopython=True)
def is_inside_sm(polygon, point):
    # From https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
    # and
    # https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1


@numba.njit(parallel=True)
def is_inside_sm_parallel(points, polygon):
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        D[i] = is_inside_sm(polygon,points[i])
    return D


def run():
    cellBoundaries_file = r"D:\rotated_dapi_map_tiles\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\cellBoundaries\cellBoundaries.tsv"
    cell_props_file = r"D:\rotated_dapi_map_tiles\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\cell_props\cell_props.tsv"
    clip_poly_file = r"D:\rotated_dapi_map_tiles\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\roi\roi_rotated.csv"

    cellBoundaries = pd.read_csv(cellBoundaries_file, sep='\t')
    cellProps = pd.read_csv(cell_props_file, sep='\t')
    roi = pd.read_csv(clip_poly_file)

    # the boundaries appear to be strings. Convert them to a list of tuples
    _cell_boundaries = [eval(d) for d in cellBoundaries.cell_boundaries]
    for i, x in enumerate(_cell_boundaries):
        _cell_boundaries[i] = [tuple(d) for d in x]
    cellBoundaries.cell_boundaries = _cell_boundaries

    # add a cell key column to cellBoundaries dataframe
    assert np.all(cellProps.cell_label == cellBoundaries.cell_label), 'cellBoundaries and cell_props are not aligned'
    cellBoundaries['cell_key'] = cellProps.cell_key
    cellBoundaries['x_cell'] = cellProps.x
    cellBoundaries['y_cell'] = cellProps.y
    cellBoundaries['cell_area'] = cellProps.cell_area
    cellBoundaries = cellBoundaries.set_index('cell_key')

    points = cellBoundaries[['x_cell', 'y_cell']].values
    poly = roi.values
    mask = is_inside_sm_parallel(points, poly)

    df = cellBoundaries[mask]
    df = df.drop(['cell_label'], axis=1)
    df = df.rename(columns={'cell_area': 'area',
                            'cell_boundaries': 'coords'})
    df = df.sort_values(['x_cell', 'y_cell'], ascending=[True, True])
    df['label'] = np.arange(df.shape[0]).astype(np.int32) + 1

    cells_df = df[['label', 'area', 'x_cell', 'y_cell']].reset_index()
    cell_boundaries_df = df[['label', 'coords']].reset_index()
    return cells_df, cell_boundaries_df


if __name__ == "__main__":
    _cells, _cellBoundaries = run()
    print('ok')
