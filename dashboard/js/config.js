function config() {
    var ini = [
        {
            name: 'default',
            roi: {"x0": 0, "x1": 20268, "y0": 0, "y1": 27352}, // Now I am confused!! (ah ok, see note below)
            imageSize: [194250, 262144],
            cellBoundaries: './cell_coords.json',
            cellData: '../demo_data/MOUSE_FULL_CORONAL/cached_results/cellData.json',
            topo: './cellBoundaries.topojson',
            // geneData: './dashboard/data/img/default/csv/Dapi_overlays.csv',
            tiles: './img/262144px/{z}/{y}/{x}.png'
        }, // 1
    ];
    var out = d3.map(ini, function (d) {
        return d.name;
    });
    return out
}

//
// ok here it is... map tiles and fov arent produced by the same image, that's why the backround image and the data plotted are misaligned (I think!).
// You padded the original image to make the longer side a multiple of 2000 whereas you just scaled it up to the correct dimension
// when you did the pyramid tiles
//