import pandas as pd
import numpy as np
from PIL import Image, ImageOps, ImageDraw
from scipy.sparse import coo_matrix, save_npz, load_npz
import os

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


def expand_z(label_image, spots, ppm):
    """
    expands the z-dimension so that all X,Y,Z have the same resolution.
    Also, it scales up the z-coords of the spots
    """
    z, h, w = label_image.shape

    k = int(ppm)
    depth = z * k

    coo = []
    for i in range(depth):
        j = divmod(i, k)[0]
        _img = label_image[j]
        coo.append(coo_matrix(_img))

    spots.z = spots.z * ppm
    return coo, spots


def run():
    width = int(5865)  # width of the original image
    height = int(7705)  # length of the original image
    z_start = int(35)
    z_end = int(46)
    ppm = 6.0121

    # # read some demo data
    spots = pd.read_csv(os.path.join(ROOT_DIR, 'pciSeq', 'data', 'B2A3', 'spots.csv'))

    label_image = np.load(os.path.join(ROOT_DIR, 'pciSeq', 'data', 'B2A3', 'B2A3_label_image.npz'))
    label_image = label_image['arr_0']  # this is already downscaled, I downsized to feed it to cellpose (ppm = 6.0121)

    v, i = np.unique(label_image, return_inverse=True)
    label_image = i.reshape(label_image.shape).astype(np.uint32)

    label_image = label_image[z_start:z_start+1]

    label_image_big = []
    for page in label_image:
        label_image_big.append(np.array(Image.fromarray(page).resize((width, height), Image.NEAREST), dtype=np.uint32))
    label_image_big = np.stack(label_image_big)

    spots_out = spots[(spots.z >= z_start) & (spots.z < z_end)]
    spots_out = spots_out.assign(z=spots_out.z-z_start)

    coo, spots_out = expand_z(label_image_big, spots_out, ppm)

    np.savez('B2A3_label_image_truncated_2.npz', coo)
    print('list of coo matrices saved at B2A3_label_image_truncated.npz')

    pd.DataFrame(spots_out).to_csv('B2A3_spots_truncated.csv', index=False)
    print('spots saved at B2A3_spots_truncated.csv')



if __name__=="__main__":
    # coo_list = np.load('my_coo_list.npz', allow_pickle=True)['arr_0']
    run()
    print('Done!')