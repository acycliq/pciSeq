function config() {
    return {
            roi: {"x0": 0, "x1": 27352, "y0": 0, "y1": 20268 },
            imageSize: [262144, 194250],
            tiles: 'https://raw.githubusercontent.com/acycliq/full_coronal_jpg_datastore/master/262144px/{z}/{y}/{x}.jpg',
            cellData: 'https://api.github.com/repos/acycliq/full_coronal_json_files/contents/data/cell_call_demo_data/mouse_full_coronal/cell_type_output/tsv/cellData_split?ref=master',
            geneData: 'https://api.github.com/repos/acycliq/full_coronal_json_files/contents/data/cell_call_demo_data/mouse_full_coronal/cell_type_output/tsv/geneData_split?ref=master',
            cellCoords: 'https://api.github.com/repos/acycliq/full_coronal_json_files/contents/data/cell_call_demo_data/mouse_full_coronal/cell_type_output/tsv/cellCoords_split?ref=master',
            class_name_separator: '.' //The delimiter in the class name string, eg if name is Astro.1, then use the dot as a separator, if Astro1 then use an empty string. It is used in a menu/control to show the class names nested under its broader name
        }
}
