# Run this by calling
#     python setup.py sdist bdist_wheel # old way to build a package
#     or
#     python -m build                   # new way to build a package

import os
from setuptools import setup, find_packages
from sys import platform


def get_static_files(root):
    out = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            out.append(os.path.join(path, name))
    return [d.strip('./pciSeq/') for d in out
            if (d.endswith('.html')
                or d.endswith('.js')
                or d.endswith('.css')
                or d.endswith('.msi')
                or d.endswith('.so')
                or d.endswith('.json')
                or d.endswith('PotreeConverter'))
            ]


install_deps = ['numpy_groupies', 'pandas', 'dask', 'scipy', 'streamlit', 'altair',
                'scikit-image', 'scikit-learn', 'tqdm', 'flask', 'fastremap',
                'numexpr', 'diplib', 'pyvips', 'natsort', 'redis', 'pytest',
                'matplotlib', 'laspy', 'jupyterlab', 'tomlkit']

version = None
with open(os.path.join('pciSeq', '_version.py'), 'r') as fid:
    for line in (line.strip() for line in fid):
        if line.startswith('__version__'):
            version = line.split('=')[1].strip().strip('\'')
            break
if version is None:
    raise RuntimeError('Could not determine version')

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pciSeq",
    version=version,
    license="BSD",
    author="Dimitris Nicoloutsopoulos",
    author_email="dimitris.nicoloutsopoulos@gmail.com",
    description="Probabilistic cell typing for spatial transcriptomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/acycliq/pciSeq",
    packages=find_packages(),
    install_requires=install_deps,
    extras_require={
        'interactive': ['matplotlib>=2.2.0', 'jupyter'],
    },
    include_package_data=True,
    package_data={'pciSeq': get_static_files(os.path.join('pciSeq', 'static'))},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)
