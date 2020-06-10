'''
splits the json file (geneData.json or cellData.json) to n number
of smaller jsons
'''

import pandas as pd
import numpy as np
import os
import glob
import config



def splitter(str_id, n):
    # x = 'geneData' # can either be geneData or cellData
    TARGET_FILE = os.path.join('./', str_id + '.json')
    OUT_DIR = os.path.join('./', str_id + '_split')

    df = pd.read_json(TARGET_FILE)
    # n = 3
    df_list = np.array_split(df, n)


    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.json')
        for f in files:
            os.remove(f)


    for i, d in enumerate(df_list):
        fname = os.path.join(OUT_DIR, str_id + '_%d.json' % i)
        d.to_json(fname,  orient='records')

    # [d.to_json() for d, i in enumerate(df_list)]
    # print(df.head())
    # df.to_json(full_path, orient='records')


N = 15
splitter('cellData', N)
splitter('geneData', N)