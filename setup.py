import setuptools
import os
from setuptools import setup

install_deps = ['numpy', 'pandas', 'sklearn',
                'numpy_groupies', 'xarray', 'numexpr',
                'diplib', 'scikit-image', 'opencv-python',
                'tqdm', 'pyvips']

version = None
with open(os.path.join('pciSeq', 'src', '_version.py'), 'r') as fid:
    for line in (line.strip() for line in fid):
        if line.startswith('__version__'):
            version = line.split('=')[1].strip().strip('\'')
            break
if version is None:
    raise RuntimeError('Could not determine version')

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pciSeq_3D",
    version=version,
    license="BSD",
    author="Dimitris Nicoloutsopoulos",
    author_email="dimitris.nicoloutsopoulos@gmail.com",
    description="Probabilistic cell typing for spatial transcriptomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/acycliq/pciSeq",
    # setup_requires=[
    #   'pytest-runner',
    #   'setuptools_scm',
    # ],
    packages=setuptools.find_packages(),
    # use_scm_version=True,
    install_requires=install_deps,
    extras_require={
        'interactive': ['matplotlib>=2.2.0', 'jupyter'],
    },
    # include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
     entry_points = {
        'console_scripts': [
          'pciSeq = pciSeq.__main__:main']
     }
)