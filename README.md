
Probabilistic Cell typing by In situ Sequencing
==============================================

A pipeline to cell type and visualise iss data. It implements the cell calling algorithm described in [[Qian, X., et al. Nature Methods (2019)]](#1) (This writeup is not finished)

![](screencast.gif)

## How to
There are three stages involved:
- [Preprocessing](#Preprocessing) 
- [Cell typing](#Cell-typing) 
- [Visualisation](#The-viewer)

The preprocessing stage prepares the data for the cell typing step which, when finished generates the flatfiles to be fed into the viewer. 
The animated screenshot above is taken from this [demo](https://acycliq.github.io/full_coronal_section/)

## Preprocessing
The main purpose is to tell whether a spot falls within a cell soma and if that happens, which cell it is.
We need as an input the ```label_image``` array of the full image and the text file describing the spots found on the image (Gene name and spot coordinates)
A ```label_image``` is an array the same size as the dapi image. Its values represent the label of the real life object (biological cell) which is shown at the corresponding location
of the dapi image. If an element of the  ```label_image``` array is zero, then the corresponding pixel is on the background. Herein the terms ```cellmap``` and ```label_image``` are
used interchangeably, with the latter being the most common in the image segmentation community.

This stage breaks up the ```label_image``` into smaller same-sized chunks/tiles/fov going from the top-left corner of the Dapi to the right and then top to bottom (see image below:)
![](preprocessing_1.jpg)

### Configuration:
To start preprocessing you need to create a ```PREPROCESSOR``` dictionary in [config.py](./config.py)  defined as follows:
```
PREPROCESSOR: 'dict'
    dictionary with keys:
    'fov_shape':
        list of length 2. The two items of the list represent the size in pixels of the x-side and y-side of the fov respectively
    'fovs_across': 
        Number of fovs along the x-axis (covering the full length of the x-side of the image) (int)
    'fovs_down':
        Number of fovs along the y-axis (covering the full length of the y-side of the image) (int)
    'spots_full':
        the path to the csv with all the spots. Headers should be  'Gene', 'x' and 'y' 
    'cellmap_full':
        the path to the label_image of the dapi
```

### Notes
 - the total length of the all fovs arranged next to another (either vertically or horizontally) can exceed the size of the corresponding side in the image
as shown in the image above. The dapi image has ```width=27352px``` and ```height=20268px``` and we set ```fovs_across=14``` and ```fovs_across=11``` totalling to 
```28000px``` and ```22000px```, assuming that each fov is square with side length ```2000px```.
- The fov doesnt have to be square. It doesnt even to be equal to the actual fov of the microscope



## Cell typing

## Visualisation
The viewer is a javascript web application running on the client side. Main tools used are 
- [Leaflet](http://leafletjs.com) on a canvas renderer
- [Leaflet](http://leafletjs.com) on a WebGL renderer (due to the sublime [PixiOverlay](https://github.com/manubb/Leaflet.PixiOverlay) class)
- [D3.js](https://d3js.org/)

It implements a [tiled web map](https://en.wikipedia.org/wiki/Tiled_web_map), a very popular technology in GIS and map services like GoogleMaps. 
The backround image is prerendered at different zoom levels and then cut into tiles. The browser, given the zoom level we want to display, fetches all the necessary tiles,
and lays them on the screen in such a manner that they compose a seamless bigger image. 

At the smallest zoom level, the whole background image fits entirely in a small single square with side length 256 pixels. Increasing zoom level by 1 doubles the map dimensions,
hence we now need 4 tiles to cover the whole image as shown in the fig below 

![](tiling_3.jpg)

The general rule is: 

<img src="https://render.githubusercontent.com/render/math?math=length (px) = width (px) = 256 * 2^{zoom} ">

The table below lists the map size for the first 10 levels
![](map_sizes.jpg)

For a full coronal slice from a mouse brain 10 zoom levels should be enough. A smaller slice, like hippocampus, would need about 6 or 7 zoom levels. 


### Configuration:
To run the viewer you need to 
- Do the map. Use  the function `tile_maker` from the [stage_image.py](./stage_image.py) module
- Pass the correct settings in [config.js](./dashboard/js/config.js). The file contains the following javascript object ![](config.js.jpg).

More specifically the properties are:
```
@property roi:  The size of the dapi image in the form {"x0": x_min, "x1": x_max, "y0": y_min, "y1": y_max}
@property imageSize: the size of the map at zoom level = 10. If for example the dapi image is 23352px-by-20268px 
                 to serve 10 zoom levels the x-side of the map has to be 262144px long. To keep proportionality
                 the y-size turns out to be 194250px
@property tiles: the path to the map tiles. This is composed by the folder path string followed by the postfix
                 '/{z}/{y}/{x}.jpg'. (assumes that map tiles are jpgs)
@property cellData: URL of the location that keeps the cellData.tsv files. It should be in the form: 
                 'https://api.github.com/repos/:owner/:repo/contents/:path?ref=branch_name',
@property geneData: URL of the location that keeps the geneData.tsv files. Same provisions as above apply here too.
@property cellCoords: URL of the location that keeps the cellCoords.tsv files. Same provisions as above apply here too.
@property class_name_separator: If the class name follows a pattern like 'superClass.subClass.subSubClass' then 
                  classes can hierecicly get nested together and collectively switched on/off by a control at the 
                  top right of the map:
```
- Configure the glyphs by editing [glyphConfig.js](./dashboard/js/glyphConfig.js). The colour is controlled by the `grlyphColor` function
at the bottom of the file and maps each `taxonomy` to a color. `taxonomy` serves a very loose grouping that brings together genes that bradly speaking tend 
to coexist in a cell type. Hence by assigning the same colour to these genes, it will be easier for the user to guess the cell types by eye balling the map and 
the patterns arising from similarly coloured spots. If you think this is not possible or useful for your own dataset, you will have to use for `taxonony` 
something appropriate for your data or even just replace the `taxonomy` values with the gene names. In this case effectively you wll be assigning a colour directly.
to a gene. Note however that 9 glyphs are avaliable. If you need more you will have to design them yourself. The glyphs shapes are declared in 
[glyphPaths.js](./dashboard/js/glyphPaths.js)
    


## References 
<a id="1">[1]</a> 
Qian, X., et al. (2019). Probabilistic cell typing enables fine mapping of closely related cell types in situ. Nat
Methods 17, 101 – 106.

## Contact
[DN](mailto:dimitris.nicoloutsopolos@gmail.com) 
