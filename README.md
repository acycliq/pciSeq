# pciSeq: Probabilistic Cell typing by In situ Sequencing
A Python package that implements the cell calling algorithm as described in [[Qian, X., et al. Nature Methods (2020)]](https://www.nature.com/articles/s41592-019-0631-4)
<p align="center">
    <img src="viewer/assets/screencast_resized.gif" alt="screenshot"/>
</p>

## Installation
```
pip install pciSeq
```
Requirement: Python >= 3.7

## Usage
You need to create two `pandas dataframes` for the spots and the single cell data and a `coo_matrix` for the label image (which in 
most cases will be the output of some image segmentation application). Then you pass them into the `pciSeq.fit()` method as follows: 
```
import pciSeq

res = pciSeq.fit(spots_df, label_image, scRNA_df)
```
See the demo below for a more detailed explanation about the arguments of  `pciSeq.fit()` and its return values.

There is also a fourth argument (optional) to override the default hyperparameter values which are initialised 
by the [config.py](https://github.com/acycliq/pciSeq/blob/master/pciSeq/config.py) module. To pass user-defined hyperparameter values, create a `dictionary` with `keys` the
hyperparameter names and `values` their new values. For example, to exclude all Npy and Vip spots you can do:

```
import pciSeq

opts = { 'exclude_genes': ['Npy', 'Vip'] }
res = pciSeq.fit(spots_df, label_image, scRNA_df, opts)
```

## Demo
You can run a pciSeq demo in google colab: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/acycliq/pciSeq/blob/master/notebooks/pciSeq.ipynb)

## Visualizer
An interactive viewer to interrogate the data runs on this [url](https://acycliq.github.io/visage/). Instructions about 
building this viewer with your own data are [here](https://github.com/acycliq/visage).

## References 
Qian, X., et al. (2020). Probabilistic cell typing enables fine mapping of closely related cell types in situ. Nat
Methods 17, 101 - 106.


