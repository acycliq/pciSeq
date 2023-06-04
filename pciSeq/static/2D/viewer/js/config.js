// NOTES:
// 1. roi is the image size in pixels. Leave x0 and y0 at zero and set x1 to the width and y1 to the height
// 2. zoomLevels, leave this at 10
// 3. layers is a dict that keeps the background image(s) of the viewer. The keys are just a label, the values
//    should point to the folder that keeps the pyramid of tiles for that particular image.
//    If you do not have that just change the link to a blind one (change the jpg extension for example).
//    The viewer should work without the dapi background though.
// 4. cellData, geneData, cellBoundaries:
//    4.1: mediaLink paths are with respect to the location of 'streaming-tsv-parser.js'
//    4.2: size is the tsv size in bytes.  Not crucial if you dont get it right, ie the full tsv will
//         still be parsed despite this being wrong. It is used by the loading page donutcharts to calc
//         how far we are.


function config() {
    return {
            roi: {"x0": 0, "x1": 7602, "y0": 0, "y1": 5471},
            zoomLevels: 10, // maximum zoom levels. Leave that at 10.
            layers: {
                'dapi': 'https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg',
                'dapi_shuffled': 'https://storage.googleapis.com/ca1-data/img/262144px/{z}/{x}/{y}.jpg',
                'empty': '',
                'world_shuffled': 'http://{s}.tile.osm.org/{z}/{y}/{x}.png',
                'world': 'http://{s}.tile.osm.org/{z}/{x}/{y}.png'
            },
            cellData: {mediaLink: '../data/cellData.tsv', size: "2179232"},
            geneData: {mediaLink: '../data/geneData.tsv', size: "9630856"},
            cellBoundaries: {mediaLink: '../data/cellBoundaries.tsv', size: "1306209"},
        }
}

