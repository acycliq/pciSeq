function config() {
    var ini = [
        {
            name: 'default',
            roi: {"x0": 0, "x1": 40000, "y0": 0, "y1": 46000},
            imageSize: [227951, 262144],
            cellBoundaries: './cellBoundaries.json',
            cellData: './cellData.json',
            topo: './cellBoundaries.topojson',
            // geneData: './dashboard/data/img/default/csv/Dapi_overlays.csv',
            tiles: 'https://raw.githubusercontent.com/acycliq/viewer_datastore/master/demo_data/human/map_tiles/262144px/{z}/{y}/{x}.png'
        }, // 1
    ];
    var out = d3.map(ini, function (d) {
        return d.name;
    });
    return out
}

