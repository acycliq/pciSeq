---
title: Command line
description: Run pciSeq from a config file with the pciseq command, override settings for a single run, and drive parameter sweeps.
---

# Command line

`pciseq` runs a dataset described by a config file. It loads the inputs and calls
[`fit`](./reference#fit), so a run is defined by a file rather than by a script.

```bash
pciseq run analysis.yaml
```

## The config file

The file has four sections: `spots`, `masks` and `scrnaseq` are required, `opts`
is optional.

```yaml
scrnaseq:
  path: scRNAseq.csv
  index_col: 0                 # column holding the gene names
  transpose: false             # true if the file is classes x genes

spots:
  path: spots.parquet
  rename: {discriminability: score, Gene: gene_name, z_stack: z_plane}
  filter: "score > 3.0 and intensity > 0.15"

masks:
  path: gcamp_seg.npz
  key: seg                     # only needed for npz holding several arrays

opts:
  voxel_size: [1, 1, 3.46]
  rTheta: 2
  save_data: true
  output_path: /path/to/output
```

Everything under `opts` goes to `fit`, so the keys are the ones on the
[configuration page](./configuration). Keys left out keep their defaults.

### scrnaseq

The cell type definitions, genes by cell types.

| key | meaning |
|-----|---------|
| `path` | csv or tsv file |
| `index_col` | column holding the gene names, by name or position. Default `0` |
| `sep` | column separator. Default `,`, use `"\t"` for tsv |
| `transpose` | `true` when the file has classes as rows |

### spots

| key | meaning |
|-----|---------|
| `path` | `.parquet`, or csv/tsv |
| `sep` | separator for text files |
| `rename` | old name to new name, for files whose columns are not what pciSeq expects |
| `filter` | a pandas query string. Several conditions join with `and` |

pciSeq needs `gene_name`, `x`, `y` and, for 3D, `z_plane`. `rename` maps the
file's column names onto these.

`filter` is applied on load, so the file on disk stays unfiltered and the
threshold stays in the config file.

### masks

| key | meaning |
|-----|---------|
| `path` | `.npy` holding the label volume, or `.npz` |
| `key` | which array to take from an npz. Required when the npz holds more than one |

A `(Z, Y, X)` volume becomes one plane per z. A 2D array is treated as a single
plane.

## Overriding settings

`--set` changes one setting without touching the file. Repeat it as needed.

```bash
pciseq run analysis.yaml --set rTheta=5 --set mrf_beta=1.0
```

Values are parsed as JSON:

```bash
--set rTheta=5                       # number
--set save_data=false                # boolean
--set 'voxel_size=[0.28,0.28,0.7]'   # list, quoted so the shell leaves the brackets alone
--set output_path=default            # not valid JSON, so kept as a string
```

A key that is not a pciSeq setting raises an error and stops the run.

## Sweeps

A sweep runs one config file in a shell loop, with one output directory per run:

```bash
for t in 2 5 15 25; do
  for b in 0 0.5 1.0; do
    pciseq run analysis.yaml \
      --set rTheta=$t --set mrf_beta=$b \
      --set output_path=out/rtheta${t}_beta${b}
  done
done
```

A change to `analysis.yaml` applies to every run of the sweep. A large 3D dataset
can take over 8 GB per run, so sweeps are best run sequentially.

## Validating a config

`--dry-run` prints the settings that would be used, as JSON, and stops without reading
any data. The `opts` are resolved as in [`fit`](./reference#fit): the defaults, then the
file's `opts`, then the `--set` values. An unknown option name or a value of the wrong
type is an error, so a typo is caught before any data is loaded.

```bash
pciseq run analysis.yaml --set rTheta=5 --dry-run
```

## Version

```bash
pciseq --version
```

The model changes between versions, so the version is part of the description of any
result. Each run also records its provenance in the output, see
[run provenance](../installation#run-provenance).

## Numbers in YAML

Numbers may be written in exponent form, `1.0e6`. PyYAML follows YAML 1.1, which
requires a signed exponent, so it reads the unsigned form as the string `"1.0e6"`; pciSeq accepts it
as a number.

## Config handling from Python

The command is a wrapper over [`fit`](./reference#fit). The config handling is
available without it:

```python
from pciSeq.cli import load_config, build_run
from pciSeq import fit

cfg = load_config("analysis.yaml")
spots, coo, scRNAseq, opts = build_run(cfg, overrides=["rTheta=5"])
cellData, geneData = fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)
```

## Beyond the config file

The config covers renaming, filtering and selecting an array from an npz. Other
preprocessing, such as a coordinate transform, merging files or stitching tiles, belongs
in a script that writes out clean inputs for the config to point at.
