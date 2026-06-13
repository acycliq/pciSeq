# Run this by calling
#     python setup.py sdist bdist_wheel # old way to build a package
#     or
#     python -m build                   # new way to build a package

import os
import subprocess
import datetime
from setuptools import setup, find_packages


def _git(args):
    """Run a git command in the source tree, return stripped stdout or ''."""
    try:
        return subprocess.check_output(
            ['git'] + args,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except Exception:
        return ''


def write_build_info():
    """Bake the current branch and commit hash into pciSeq/_build_info.py.

    This file is gitignored and regenerated on every build, so the values
    end up in the wheel and survive pip install (where .git is gone).
    """
    commit = _git(['rev-parse', '--short', 'HEAD'])
    branch = _git(['rev-parse', '--abbrev-ref', 'HEAD'])
    build_date = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    path = os.path.join('pciSeq', '_build_info.py')
    with open(path, 'w') as f:
        f.write(f'__commit__ = {commit!r}\n')
        f.write(f'__branch__ = {branch!r}\n')
        f.write(f'__build_date__ = {build_date!r}\n')


write_build_info()


def get_static_files(root):
    out = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            out.append(os.path.join(path, name))
    return [
        d.strip("./pciSeq/")
        for d in out
        if (
            d.endswith(".html")
            or d.endswith(".js")
            or d.endswith(".css")
            or d.endswith(".so")
            or d.endswith(".json")
            or d.endswith("PotreeConverter")
        )
    ]


install_deps = [
    "numpy_groupies",
    "pandas",
    "dask",
    "scipy",
    "scikit-learn",
    "tqdm",
    "flask",
    "flask-socketio",
    "fastremap",
    "diplib",
    "pyvips[binary]",
    "natsort",
    "matplotlib",
    "colorlog",
    "shapely",
    "alphashape",
    "opt_einsum",
    "plotly",
    "numba",
    "pyarrow",
]


def get_version():
    """Get version from _version.py and append git commit hash if available."""
    version = None
    with open(os.path.join("pciSeq", "_version.py"), "r") as fid:
        for line in (line.strip() for line in fid):
            if line.startswith("__version__"):
                version = line.split("=")[1].strip().strip("'\"")  # Strip both ' and "
                break

    if version is None:
        raise RuntimeError("Could not determine version")

    # Try to append git commit hash
    # NOTE: Disabled because +g{commit} format causes pip install failures in CI
    # The version string must be PEP 440 compliant
    # try:
    #     import subprocess
    #     commit = subprocess.check_output(
    #         ['git', 'rev-parse', '--short', 'HEAD'],
    #         stderr=subprocess.DEVNULL,
    #         text=True
    #     ).strip()
    #     version = f"{version}+g{commit}"
    # except Exception:
    #     # No git or not in a git repo - use version as is
    #     pass

    return version


version = get_version()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pciSeq_3d",
    version=version,
    license="BSD",
    author="Dimitris Nicoloutsopoulos",
    author_email="dimitris.nicoloutsopoulos@gmail.com",
    description="Probabilistic cell typing for spatial transcriptomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/acycliq/pciSeq_3d",
    packages=find_packages(),
    install_requires=install_deps,
    include_package_data=True,
    package_data={
        "pciSeq": get_static_files(os.path.join("pciSeq", "static"))
        + get_static_files(os.path.join("pciSeq", "src", "realtime_viewer"))
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)
