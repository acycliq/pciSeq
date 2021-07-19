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


def stage_vizgen():
    cellBoundaries_file = r"D:\rotated_dapi_map_tiles\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\cellBoundaries\cellBoundaries.tsv"
    cell_props_file = r"D:\rotated_dapi_map_tiles\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\cell_props\cell_props.tsv"
    spots_file = r"D:\Home\Dimitris\OneDrive - University College London\dev\Python\pciSeq\pciSeq\data\vizgen\merfish\labelled_spots.tsv"
    clip_poly_file = r"D:\rotated_dapi_map_tiles\MsBrain_Eg1_VS6_JH_V6_05-02-2021\region_0\roi\roi_rotated.csv"

    cellBoundaries = pd.read_csv(cellBoundaries_file, sep='\t')
    cellProps = pd.read_csv(cell_props_file, sep='\t')
    roi = pd.read_csv(clip_poly_file)
    spots_df = pd.read_csv(spots_file, sep='\t')

    # the boundaries appear to be strings. Convert them to a list of tuples
    _cell_boundaries = [eval(d) for d in cellBoundaries.cell_boundaries]
    for i, x in enumerate(_cell_boundaries):
        _cell_boundaries[i] = [tuple(d) for d in x]
    cellBoundaries.cell_boundaries = _cell_boundaries

    # add a cell key column to cellBoundaries dataframe
    assert np.all(cellProps.cell_label == cellBoundaries.cell_label), 'cellBoundaries and cell_props are not aligned'
    cellBoundaries['cell_key'] = cellProps.cell_key
    cellBoundaries['x'] = cellProps.x
    cellBoundaries['y'] = cellProps.y
    cellBoundaries['cell_area'] = cellProps.cell_area
    cellBoundaries = cellBoundaries.set_index('cell_key')

    points = cellBoundaries[['x', 'y']].values
    poly = roi.values
    mask = is_inside_sm_parallel(points, poly)

    df = cellBoundaries[mask]
    df = df.drop(['cell_label'], axis=1)
    df = df.rename(columns={'cell_area': 'area',
                            'cell_boundaries': 'coords'})

    # for some reason two or more cells can have the same cell_key. I havent
    # looked at it but for now I will remove the duplicates and keep the one
    # with the largest area
    df = df.sort_values(by=['cell_key', 'area'], ascending=[False, False])
    df = df[~df.index.duplicated(keep='first')]

    # after this data cleaning, sort by the centroid
    df = df.sort_values(['x', 'y'], ascending=[True, True])
    df['label'] = np.arange(df.shape[0]).astype(np.int32) + 1

    cells_df = df[['label', 'area', 'x', 'y']].reset_index().drop(['cell_key'], axis=1)
    cell_boundaries_df = df[['label', 'coords']].reset_index()

    spots_df = spots_df.drop(['neighbour', 'neighbour_array', 'neighbour_prob'], axis=1)
    spots_df = spots_df.rename(columns={'x': 'x_global',
                                        'y': 'y_global',
                                        'Gene': 'target'})
    spots_df.inside_cell_key = spots_df.inside_cell_key.fillna(0)

    spots_df = spots_df.merge(df, how='left', left_on='inside_cell_key', right_index=True)
    spots_df = spots_df[['x_global', 'y_global', 'label', 'target', 'x', 'y']]
    spots_df.label = spots_df.label.fillna(0)
    spots_df = spots_df.rename(columns={'x': 'x_cell',
                                        'y': 'y_cell'})
    spots_df = spots_df.astype({'x_global': np.int32,
                                'y_global': np.int32,
                                'label': np.int32})



    return cells_df, cell_boundaries_df, spots_df


if __name__ == "__main__":
    _cells, _cellBoundaries, spots_df = stage_vizgen()
    print('ok')
