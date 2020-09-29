import pandas as pd
import numpy as np
import cv2
from scipy import ndimage
from skimage import morphology
from skimage import exposure, feature
from skimage.measure import regionprops
from skimage.segmentation import watershed
from scipy.ndimage import label, generate_binary_structure
from scipy.sparse import coo_matrix
from src.preprocess.imimposemin import imimposemin
from skimage.color import label2rgb
from matplotlib import pyplot as plt


# set the parameters
DapiThresh = 90
DapiMinSize = 5
DapiMinSep = 7
DapiMargin = 10
MinCellArea = 200

def imadjust(img):
    upper = np.percentile(img, 99)
    lower = np.percentile(img, 1)
    out = (img - lower) * (255 / (upper - lower))
    np.clip(out, 0, 255, out) # in-place clipping
    return np.around(out)

def disk(n):
    struct = np.zeros((2 * n + 1, 2 * n + 1))
    x, y = np.indices((2 * n + 1, 2 * n + 1))
    mask = (x - n)**2 + (y - n)**2 <= n**2
    struct[mask] = 1
    return struct.astype(np.uint8)

def imregionalmax(image, ksize=3):
  """Similar to matlab's imregionalmax"""
  filterkernel = np.ones((ksize, ksize)) # 8-connectivity
  reg_max_loc = feature.peak_local_max(image,
                               footprint=filterkernel, indices=False,
                               exclude_border=0)
  return reg_max_loc.astype(np.uint8)

# selem = morphology.disk(20)
# morphology.erosion(image, selem)

Dapi = cv2.imread(r"C:\Users\skgtdni\Downloads\block_18.tif", cv2.IMREAD_GRAYSCALE)

Dapi = imadjust(Dapi)
ThresVal = np.percentile(Dapi[Dapi>0], DapiThresh)
kernel = disk(2)
# image = cv2.erode(Dapi>ThresVal), kernel)
bwDapi = ndimage.binary_erosion(Dapi>ThresVal, structure=disk(2)).astype(int)
dist = ndimage.distance_transform_edt(bwDapi)
dist0 = dist
dist0[dist<DapiMinSize]=0

selem = morphology.disk(DapiMinSep)
ddist = morphology.dilation(dist0, selem)

impim = imimposemin(-dist0, imregionalmax(ddist))

bwDapi0 = bwDapi
bwDapi0[watershed(impim) == 0] = 0

# % assign all pixels a label
s = generate_binary_structure(2,2)
labels, num_Features = label(bwDapi0, structure=s)
d, _idx = ndimage.distance_transform_edt(bwDapi == 0, return_indices=True)
idx = np.ravel_multi_index(_idx, bwDapi.shape)


# % now expand the regions by a margin
CellMap0 = np.zeros(Dapi.shape).astype(np.uint32)
Expansions = d < DapiMargin
# CellMap0[Expansions] = labels[idx[Expansions]];
CellMap0 = np.take(labels.flatten(), idx) * Expansions

rProps0 = regionprops(CellMap0)

# BigEnough = np.array([x.area > 200 for x in rProps0 ])
# NewNumber = np.zeros(len(rProps0))
# NewNumber[~BigEnough] = 0
# NewNumber[BigEnough] = np.arange(1, 1+BigEnough.sum())
# CellMap = CellMap0
# CellMap[CellMap0 > 0] = NewNumber[CellMap0[CellMap0 > 0]]
coo = coo_matrix(CellMap0)
coo_data = coo.data.copy()
to_remove = np.array([x.label for x in rProps0 if x.area <= MinCellArea])
coo_data[np.in1d(coo_data, to_remove)] = 0
_, res = np.unique(coo_data, return_inverse=True)
coo.data = res

my_image = cv2.imread(r"C:\Users\skgtdni\Downloads\block_18.tif", cv2.IMREAD_GRAYSCALE)
# overlay = label2rgb(coo.toarray(), image=my_image, bg_label=0)
overlay = label2rgb(coo.toarray(), bg_label=0, bg_color=(1,1,1))
# plt.imshow(overlay[:, :, 1], cmap='gray', interpolation='none')
my_dpi = 72
fig, ax = plt.subplots(figsize=(6000/my_dpi, 6000/my_dpi), dpi=my_dpi)
plt.imshow(overlay)
ax.set_axis_off()
plt.tight_layout()
plt.show()

fig.savefig('block_18_segmented.tif', dpi=200)



print(coo.data.max())
print('done')

# https://stackoverflow.com/questions/5260232/matlab-octave-bwdist-in-python-or-c

# v = np.ones_like(i)
# mat = coo_matrix((v, (i, j)), shape=(n, n))