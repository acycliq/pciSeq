from pciSeq.src.core.logger import logger_setup
from pciSeq import fit
from scipy.sparse import load_npz
import pandas as pd


# set up the logger
logger_setup()

# read some demo data
spots = pd.read_csv('test_spots.csv')

coo = load_npz('test_label_image_coo.npz')

scRNAseq = pd.read_csv('test_scRNAseq.csv', header=None, dtype=object)
scRNAseq.columns = scRNAseq.iloc[0]
scRNAseq = scRNAseq.iloc[1:, :]
scRNAseq = scRNAseq.set_index('gene_name')
scRNAseq = scRNAseq.applymap(lambda x: eval(x))

# main task
# _opts = {'max_iter': 10}
_opts = {'save_data': False, 'launch_viewer': False, 'launch_diagnostics': False}
cellData, geneData = fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=_opts)
