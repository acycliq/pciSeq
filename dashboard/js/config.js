function config() {
    var ini = [
        {
            name: 'default',
            roi: {"x0": 0, "x1": 27352, "y0": 0, "y1": 20268 },
            imageSize: [262144, 194250],
            // cellData: './data/cell_call_demo_data/mouse_full_coronal/cell_type_output/cellData.json',
            // tiles: 'https://raw.githubusercontent.com/acycliq/full_coronal_jpg_datastore/master/262144px/{z}/{y}/{x}.jpg',
            tiles: './dashboard/img_landscape/262144px/{z}/{y}/{x}.jpg',
            // spot_json: function(d){ return "https://raw.githubusercontent.com/acycliq/full_coronal_json_files/master/data/cell_call_demo_data/" + 'mouse_full_coronal/' + 'cell_type_output/' + 'geneData_landscape_split/' + 'geneData_landscape_' + d + '.json'},
            // cell_json: function(d){ return "https://raw.githubusercontent.com/acycliq/full_coronal_json_files/master/data/cell_call_demo_data/" + 'mouse_full_coronal/' + 'cell_type_output/' + 'cellData_landscape_split/' + 'cellData_landscape_' + d + '.json'},
            gene_tsv: function(d){ return './dashboard/img/tsv/geneData_split/' + 'geneData_' + d + '.tsv'},
            cell_tsv: function(d){ return './dashboard/img/tsv/cellData_split/' + 'cellData_' + d + '.tsv'},
            cellCoords_tsv: function(d){ return './dashboard/img/tsv/cellCoords_split/' + 'cellCoords_' + d + '.tsv'},
            num_tsvs: 4, // number of tsv splits
            class_name_separator: '' //The delimiter in the class name string, eg if name is Astro.1, then use the dot as a separator, if Astro1 then use an empty string. It is used in a menu/control to show the class names nested under its broader name
        }, // 1
    ];
    var out = d3.map(ini, function (d) {
        return d.name;
    });
    return out
}
