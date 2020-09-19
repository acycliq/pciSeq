
Probabilistic Cell typing by In situ Sequencing
==============================================

A pipeline to cell type and visualise iss data. (This writeup is not finished)

![](screencast.gif)

## HOW TO
There are three stages involved:
- Preprocessing 
- Cell typing
- Visualisation

The preprocessing stage prepares the data for the cell typing step which, when finished generates the flatfiles to be fed into the viewer.

## Preprocessing
The main purpose is to tell whether a spot falls within a cell soma and if that happens, which cell it is.
We need as an input the ```label_image``` array of the full image and the text file describing the spots found on the image (Gene name and spot coordinates)
A ```label_image``` is an array the same size as the dapi image. Its values represent the label of the real life object (biological cell) which is shown at the corresponding location
of the dapi image. If an element of the  ```label_image``` array is zero, then the corresponding pixel is on the background. Herein the terms ```cellmap``` and ```label_image``` are
used interchangeably, with the latter being the most common in the image segmentation community.
This stage breaks up the ```label_image``` into smaller same-sized chunks/tiles/fov going from the top-left corner of the Dapi to the right and then top to bottom (see image below:)
![](preprocessing_1.jpg)

```javascript
Settings dictionary with keys:
    'fov_shape':
    'fovs_across'
    'fovs_down': 
    'spots_full': 
    'cellmap_full':
```

## Cell typing

## The viewer 
The viewer is a javascript web application running on the client side. Main technologies used are 
- [Leaflet](http://leafletjs.com) on a canvas renderer
- [Leaflet](http://leafletjs.com) on a WebGL renderer (due to the sublime [PixiOverlay](https://github.com/manubb/Leaflet.PixiOverlay) class)
- [D3.js](https://d3js.org/)

## Contact
[DN](mailto:dimitris.nicoloutsopolos@gmail.com) 

https://acycliq.github.io/full_coronal_section/

