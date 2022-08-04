# Visage
An interactive web-based viewer to visualise 2D spatial transcriptomics data. A demo using 
CA1 data from [Qian, X., et al. Nature Methods (2020)](https://www.nature.com/articles/s41592-019-0631-4) runs
 [here](https://acycliq.github.io/visage/)
<p align="center">
    <img src="viewer/assets/screencast_resized.gif" alt="Your image title"/>
</p>

## Instructions
The easiest way is to run `pciSeq.fit()` passing a dict with `save_data` set to `True` as the options arg:
```python
import pciSeq

opts = { 'save_data': True }
res = pciSeq.fit(spots_df, label_image, scRNA_df, opts)
```
That should save three tsv files (`geneData.tsv`, `cellData.tsv`, `cellBoundaries.tsv`) in your temp directory

#### Step 1
Download/clone the main branch from this repo and replace the three tsv files under `visage/viewer/data/` with your own

#### Step 2
In the configuration file `visage/viewer/js/config.js` you need to edit
* `roi`: This describes the dimensions of your image. Change `x1` and `y1` only and set them to the width and height 
respectively (in pixels) of your image. 
Do not change `x0`and `y0` and leave them as zeros.ll data has been m
* `zoomLevels`: Leave that to 10
* `tiles`: this is related to the background image of the viewer. However it is not necessary, it is optional. The viewer 
works without a background image. I will explain how to do that step at a later stage, hence for now just use a blind link there; 
For example, change the extension and use something like `https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg_ZZZ`
 * `cellData`: If you have followed step 2 above you will not have to change the `mediaLink` value. Just update size to the size in bytes of 
 your own `cellData.tsv`
 * `geneData`: Same as above
 * `cellBoundaries`: Same as above
 
#### Step 3
Set your own colours to the cell classes. To do this just edit `visage/viewer/js/classColors.js`. Make sure that 
all the cell class from your single cell data have been associated with a color. If you miss one, then the 
viewer will not run. 
Note also that the every cell class is bucketed under a wider type, called `IdentifiedType`. Then each `IdentifiedType` has been 
assigned a color. Hence all cell classes with the same `IdentifiedType` will have the same color.

#### Step 4
The spot colors of the spots are set inside the script `visage/viewer/js/glyphConfig.js`. This is also where you can set the shape of the glyph.
As you zoom-in and you get closer to the scene the dots will take some shape and the combination shape-color uniquely identifies a gene. 
The glyph shapes are defined in `visage/viewer/js/glyphPaths.js`. No need to edit `glyphPaths.js` unless you want to design a new shape. 
You only have to edit `glyphConfig.js` so that it describes the color scheme of your preference for your data. Just make sure every single 
gene of your gene panel is included `glyphConfig.js`; If even one is missing, the viewer will not load.

#### Step 5
Do the background image. This is not one single image but a mosaic of several small tiles, conceptually similar to the way google maps works for example. 
I will cover that step some other day, the viewer works without a background image. However if you want to add that send me an email and I will get back to you.

#### Comments
* To load the viewer you need to run `index.html' located under the project root directory. Most modern IDEs like pyCharm or VSCode will run
a built-in http server and load the html file without any hassle. Otherwise you need to serve the root dir by cd-ing into it and running the python 3 command

   `python -m http.server 8000`
   
   Then type `localhost:8000` into your browser and the viewer should load.

* Depending on your data, the dot size might sometimes turn to be too large. Let me know if you have any problems to tell you how to fine-tune the marker size.

* There is a collapsible control at the top-right of the viewer that shows a nested tree of the classes/subclasses. It is assumed that class names 
are dot separated, for example `Astro.1` or `Sst.Npy.Cort`. Nesting will not work if they are not dot-separated (although there is code to specify the separator but hasnt been tested. 
Let me know if you fall in this case). The control will just show a list of all the cell classes (un-nested)
but interactivity will still work and the user should be able to hide/show classes by toggling the checkboxes.