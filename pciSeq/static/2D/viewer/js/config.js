// NOTES:
// 1. roi is the image size in pixels. Leave x0 and y0 at zero and set x1 to the width and y1 to the height
// 2. zoomLevels, leave this at 10
// 3. tiles should point to the folder that keeps your pyramid of tiles. If you do not have that just
//    change the link to a blind one (change the jpg extension for example). The viewer should work
//    without the dapi background though
// 4. cellData, geneData, cellBoundaries:
//    4.1: mediaLink paths are with respect to the location of 'streaming-tsv-parser.js'
//    4.2: size is the tsv size in bytes. I use os.path.getsize() to get it. Not crucial if you dont get it right,
//         ie the full tsv will still be parsed despite this being wrong. It is used by the loading page
//         donut charts to calc how far we are.

function config() {
    return {
            roi: {"x0": 0, "x1": 7602, "y0": 0, "y1": 5471},
            zoomLevels: 10, // maximum zoom levels. Leave that at 10.
            tiles: 'https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg',
            cellData: {mediaLink: '../data/cellData.tsv', size: "2179232"},
            geneData: {mediaLink: '../data/geneData.tsv', size: "9630856"},
            cellBoundaries: {mediaLink: '../data/cellBoundaries.tsv', size: "1306209"},
        }
}

