function myUtils() {
    function poly_collection(data, t) {
        console.log('Doing geojson object for polygons');
        var dots = {
            type: "FeatureCollection",
            features: []
        };
        for (var i = 0; i < data.length; ++i) {
            var c = data[i].coords,
                temp = [];
            // if (data[i].cell_id === 737) {
            //     console.log('stop')
            // }
            if (c) {
                // That needs some attention. It leaves an open bug
                // if c is null (ie a cell doesnt have boundary coords for some reason) then
                // the geometry g which will be attached to data[i] will be the previous point's (ie data[i-1]) geometry,
                // I am turning a blind eye cause I dont have too many such cells (in fact all cells should have boundaries)
                // but it is a bug!
                for (var j = 0; j < c.length; ++j) {
                    var x = c[j][0],
                        y = c[j][1];
                    var lp = t.transform(L.point([x, y]));
                    temp.push([lp.x, lp.y])
                }
                var g = {
                    "type": "Polygon",
                    "coordinates": [temp]
                };
            }

            var target_cell = get_cell(+data[i].cell_id, i);
            //create feature properties
            var p = {
                // "fov": getFov(+target_cell.X, +target_cell.Y), // <-- Do not use that if you can. Browser will load a bit faster if removed
                "id": i,
                "cell_id": +data[i].cell_id,
                "centroid": get_centroid(target_cell),
                "X": +target_cell.X,
                "Y": +target_cell.Y,
                "Genenames": target_cell.Genenames,
                "CellGeneCount": target_cell.CellGeneCount,
                "ClassName": target_cell.ClassName,
                "topClass": target_cell.topClass,
                "Prob": target_cell.Prob,
                "agg": target_cell.agg,
            };

            //create features with proper geojson structure
            dots.features.push({
                "geometry": g,
                "type": "Feature",
                "properties": p
            });
        }
        console.log('geojson done');
        return dots;
    }

    function get_cell(cell_num, i) {
        // that can be done better i think!
        if (cellData[i].Cell_Num === cell_num){
            return cellData[i]
        }
        else{
            console.log('"cellData" and "data" arrays arent aligned');
            return null
        }
        //
        // That is a safer way to do the same but takes longer
        // return cellData.filter(d => d.Cell_Num === cell_num)[0];
    }

    function get_centroid(target_cell) {
        // var target_cell = cellData.filter(d => d.Cell_Num === cell_num)[0];
        return [target_cell.X, target_cell.Y]
    }

    function getFov(x, y) {
        var turfPoint = turf.point([x, y]);

        var p = block_boundaries.features.filter(d => turf.booleanPointInPolygon(turfPoint, turf.polygon(d.geometry.coordinates)));
        if (p) {
            return p[0].properties.id
        } else {
            return null
        }
    }

    /**
     * Converts a hexadecimal string to a hexadecimal color number.
     *
     * @example
     * PIXI.utils.string2hex("#ffffff"); // returns 0xffffff
     * @memberof PIXI.utils
     * @function string2hex
     * @param {string} The string color (e.g., `"#ffffff"`)
     * @return {number} Number in hexadecimal.
     */
    function string2hex(string) {
        if (typeof string === 'string' && string[0] === '#') {
            string = string.substr(1);
        }
        return parseInt(string, 16);
    }

    function stripper(d, k) {
        for (i = 0; i < k; ++i) {
            if (d.lastIndexOf('.') > 0) {
                d = d.substring(0, d.lastIndexOf('.'))
            }
        }
        return d
    }

    // find the position of the i-th occurrence of substring m in string str
    function getPosition(str, m, i) { return str.split(m, i).join(m).length; }

    // fw_stripper('this.is.a.test', 2) = 'this.is'
    // fw_stripper('this.is.a.test', 3) = 'this.is.a'
    function fw_stripper(d, k) {
        var out = d.substring(0, getPosition(d, '.', k));
        return out
    }


    var res = {};
    res.poly_collection = poly_collection;
    res.string2hex = string2hex;
    res.stripper = stripper;
    res.fw_stripper = fw_stripper

    return res
}

