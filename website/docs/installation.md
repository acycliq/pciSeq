# Installation

Install from the GitHub repository:

```bash
pip install git+https://github.com/acycliq/pciSeq_3d.git@dev_3d
```

This installs pciSeq, including the live viewer and its chat.

There is one optional extra, `mcp`. It adds `pciseq-mcp`, the program through which an AI
assistant such as Claude Code reads the files a completed run saved (see
[MCP server](./api/mcp-server.md)). It is required only for that means of access; the other
three do not use it: an AI assistant connected to pciSeq Viewer, the chat panel in
pciSeq Viewer, and the chat of the live viewer.

```bash
pip install "pciSeq_3d[mcp] @ git+https://github.com/acycliq/pciSeq_3d.git@dev_3d"
```

For development, clone it and install in editable mode so edits to the source are picked up
without reinstalling. The extra is given after the dot:

```bash
git clone https://github.com/acycliq/pciSeq_3d.git
cd pciSeq_3d
pip install -e .           # pciSeq
pip install -e ".[mcp]"    # adds pciseq-mcp
```

Python 3.10 or later. pip installs the dependencies. The image tiling functions require
libvips, which is supplied by the `pyvips[binary]` wheel; on a platform with no such wheel the
rest of the package still imports and only `tile_maker` and `stage_image` raise an ImportError, with a
message pointing at the libvips install page. `read_tiles` does not need libvips.

## Versions

```python
import pciSeq
pciSeq.__version__     # '0.0.66.dev0'
pciSeq.__commit__      # short git hash, eg '6db6dc8'
pciSeq.__branch__      # 'dev_3d'
pciSeq.__commit_date__  # when that commit was made, eg '2026-09-25T15:02:11+01:00'
```

The version string lives in `pciSeq/_version.py`. `__commit__`, `__branch__` and
`__commit_date__` are read from git when the package is imported from a checkout,
otherwise from the values recorded when it was built. On the command line, `pciseq --version` prints the version string.

The model changes between commits. Two results are directly comparable only when the
commit matches, or when the changes in between are known. Every run records the commit; see [Run provenance](#run-provenance).

## Run provenance

Every run records the code and the environment that produced it, so a result can be
traced back to them:

```python
varBayes.metadata
# {'version':          '0.0.66.dev0',
#  'branch':           'dev_3d',
#  'commit':           '6db6dc8',
#  'commit_date':      '2026-09-04T11:30:02+01:00',
#  'created_at':       '2026-09-04T12:47:51Z',
#  'python_version':   '3.10.18',
#  'os':               'Linux 6.8.0-138-generic',
#  'package_versions': {'numpy': '2.2.6', 'scipy': '1.15.3',
#                       'pandas': '2.3.3', 'numba': '0.62.1'}}
```

No hostname or user name is recorded.

`varBayes` is the fitted model returned by [`cell_type`](api/reference.md#cell-type), or
loaded from the pickle when `save_data` is on. The same dictionary is written into the saved
output in two more places, so it survives without the live object:

- `diagnostics.db`, under the `pciSeq_provenance` key of the `metadata` table:

  ```python
  import sqlite3, json

  con = sqlite3.connect('<output_path>/pciSeq/data/viewer_data/diagnostics/diagnostics.db')
  row = con.execute(
      "select value from metadata where key = 'pciSeq_provenance'").fetchone()
  print(json.loads(row[0]))
  ```

- the [SpatialData store](api/spatialdata-store.md#provenance), as `sdata.attrs['pciseq']`,
  together with the resolved run config, without `label_map`.

`commit` is the field to quote when comparing two runs.
