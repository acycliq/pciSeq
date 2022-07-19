import pandas as pd
import numpy as np
from PIL import Image, ImageOps, ImageDraw
import skimage.io
from scipy.sparse import coo_matrix
import numpy as np
import os
import cv2

ppm = 6.01  # resolution, pixels per micron
width = 5865  # width of the original image
height = 7705  # length of the original image

# get the downscaled image
label_image = np.load(r"D:\Home\Dimitris\OneDrive - University College London\dev\Python\pciSeq\pciSeq\data\B2A3\B2A3_label_image.npz")
label_image = label_image['arr_0']

z, h, w = label_image.shape

k = int(ppm)
depth = z*k
temp_img = np.zeros([depth, h, w])

out = []
for i in range(temp_img.shape[0]):
    j = divmod(i, k)[0]
    _img = np.array(Image.fromarray(label_image[j]).resize((width, height), Image.NEAREST), dtype=np.uint32)
    out.append(coo_matrix(_img))
    # temp_img[i] = label_image[j]

print('ok')

