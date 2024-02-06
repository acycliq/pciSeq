// NOTES:
// 1. paths in 'cellData', 'geneData' and 'cellBoundaries' are with respect to the location of
//    'streaming-tsv-parser.js'
// 2. size is the tsv size in bytes. I use os.path.getsize() to get it. Not crucial if you
//    don't get it right, ie the full tsv will still be parsed despite this being wrong. It
//    is used by the loading page piecharts to calc how far we are.
// 3. roi is the image size in pixels. Leave x0 and y0 at zero and set x1 to the width and y1 to the height.
// 4. layers is a dict. Each key/value pair contains the string (the name) of the background image and the
//    location of the folder that the corresponding pyramid of tiles. If the tiles are stored locally, they
//    should be kept in a folder which is served, for example next to the tsv flatfiles. The path should be
//    in relation to the location of the index.html If you do not have a pyramid of tiles just
//    change the link to a blind one (change the jpg extension for example or just use an empty string).
//    The viewer should work without the dapi background though.
//    If the dict has more than one entries then a small control with radio button will appear at the top
//    right of the viewer to switch between different background images.
// 5. maxZoom: maximum zoom levels. In most cases a value of 8 if good enough. If you have a big image, like
//    full coronal section for example then a value of 10 would make sense. Note that this should be typically
//    inline with the zoom level you used when you did
//    the pyramid of tiles. No harm is it less. If it is greater, then for these extra zoom levels there will
//    be no background image.
// 6. spotSize: Scalar. Use this to adjust the screen-size of your spots before they morph into glyphs.
function config() {
    return {
        "cellData": { "mediaLink": "../../data/cellData.tsv", "size": "2180603"},
        "geneData": { "mediaLink": "../../data/geneData.tsv", "size": "9630820"},
        "cellBoundaries": {"mediaLink": "../../data/cellBoundaries.tsv", "size": "1306209"},
        "roi": {"x0": 0, "x1": 7602, "y0": 0, "y1": 5471},
        "maxZoom": 8,
        "layers": {
            // "empty": "",
            "dapi": "https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg"
        },
        "spotSize": 1/16
    }
}
