import os
import platform

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
print(ROOT_DIR)


PREPROCESSOR = {
    'MATLAB_SPOTS': os.path.join(ROOT_DIR, 'data', 'from_Matlab', 'split', 'spots'),
    'FOV_ROOT': os.path.join(ROOT_DIR, 'data', 'fov'),
    'fov_size': 2000,   # implies that each fov is square with side length 2000px
    'FOVS_ACROSS': 11,
    'FOVS_DOWN': 14,
    # 'fovs_across': 11,  # This is not an input, it is calculated in '''split_image'''. I put it here for convenience
    # 'fovs_down': 14,    # This is not an input, it is calculated in '''split_image'''. I put it here for convenience
}
