"""
compares the data calling results with these obtained when scRNAseq data were available.
"""
import pandas as pd
import numpy as np
import json
import os

cellData_path_2d = r'https://www.googleapis.com/download/storage/v1/b/ca1-data/o/cellData%2FcellData.tsv?generation=1628591860441109&alt=media'
celldata_path_3d = r"https://www.googleapis.com/download/storage/v1/b/ca1_3d/o/cellData.tsv?generation=1645141005940382&alt=media"

def run():
    cellData_2d = pd.read_csv(cellData_path_2d, sep='\t')
    cellData_3d = pd.read_csv(celldata_path_3d, sep='\t')

    # # get the gene panel (sorted)
    # genes = np.unique(ca1_geneData.Gene)
    # class_names = [d+'_class' for d in genes]
    # class_names.append('Zero')
    # m = cellData.shape[0]
    # n = len(class_names)
    # out = pd.DataFrame(np.zeros([m, n]), columns=class_names)


    actual_list = []
    for i, row in cellData_2d.iterrows():
        # 1. Find the most likely class as this is given by the ca1 data with proper scRNAseq
        arr = json.loads(row.Prob)
        argmax = np.argmax(arr)
        ca1_class = eval(row.ClassName)
        actual = ca1_class[argmax]
        actual_list.append(actual)

        # 2 For the same cell find the cell call results under unknown scRNAseq data
        # Do some sanity checking first
        assert row['X'] == cellData_3d.X_0[i]
        assert row['Y'] == cellData_3d.Y_0[i]

        cols = eval(cellData_2d.ClassName[i])
        vals = eval(cellData_2d.Prob[i])
        out.iloc[i][cols] = vals

    # make an extra column with the actual classes
    out['actual'] = actual_list
    # set it as the index of the df
    out = out.set_index('actual')
    # group by actual class
    res = out.groupby(out.index.values).agg('mean')
    # save the confusion matrix to csv
    res.to_csv('confusion_matrix.csv')
    print('Done')


if __name__ == "__main__":
    run()