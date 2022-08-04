import numpy as np
import pandas as pd


# def truncate_data(spots, label_image, z_min, z_max):
#     spots = truncate_spots(spots, z_min, z_max)
#     label_image = trunctate_zstack(label_image, z_min, z_max)
#     return spots, label_image
#
#
# def trunctate_zstack(masks, z_min, z_max):
#     return masks[z_min: z_max+1]
#
#
# def truncate_spots(spots, zmin, zmax):
#     spots_min = spots[(spots.z <= zmax) & (spots.z >= zmin)]
#     spots_min.assign(z=spots_min.z - zmin)
#     # out = spots_min.z - zmin
#     return spots_min


if __name__=="__main__":
    zstack_min = 18
    zstack_max = 43
    spots = pd.read_csv(r'E:\data\Anne\220308 50umCF seq atto425 DY520XL MS002\spots_yxz.csv')
    masks = np.load(r'E:\data\Anne\220308 50umCF seq atto425 DY520XL MS002\masks_2D_stiched_fullsize.npz')['arr_0']

    truncate_data(spots, masks, zstack_min, zstack_max)
    print('All done')