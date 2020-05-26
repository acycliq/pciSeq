function config() {
    var ini = [
        {
            name: 'default',
            roi: {"x0": 0, "x1": 20268, "y0": 0, "y1": 27352}, // Now I am confused!! (ah ok, see note below)
            imageSize: [194250, 262144],
            cellBoundaries: './dashboard/cell_coords.json',
            cellData: './data/cell_call_demo_data/mouse_full_coronal/cell_type_output/cellData.json',
            tiles: 'https://raw.githubusercontent.com/acycliq/full_coronal_datastore/master/{z}/{y}/{x}.png',
            spot_json: function(d){ return "./data/fov/" + 'fov_' + d + '/cell_type_out/fov_' + d + '_Dapi_overlays.json'},
            cell_json: function(d){ return "./data/fov/" + 'fov_' + d + '/cell_type_out/fov_' + d + '_iss.json'},
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
