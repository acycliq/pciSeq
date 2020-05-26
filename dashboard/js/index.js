//
// General remark. Inconsistent notation all over the place. Cell_Num and cell_id, id are the same and actually
// should be removed and replaced by label (or cell_label whatever name you decide to use)
//

// Variables in the global scope
var cellBoundaries,
    cellData,
    all_geneData = [],
    dapiConfig,
    polygonsOverlay,
    cellPolygons,
    cellSpritesLayer,
    configSettings,
    cellSpitesLayer,
    cellPolyLayer,
    glyphToNeighbourLayer,
    nearbyCellLayer,
    dotLayer,
    myDots,
    cellWatch, //keeps the id of the cell currently drawn on the map
    map,
    container_array = [],
    pixiContainer,
    pixiRenderer,
    // layerControl, //you sure this has to be global?
    geneLayers,
    geneOverlays,
    legendWindow,
    legend_added = false; //used to make sure the listener is attached only once
    pinnedControls = false;


localStorage.clear();


// this will be watching if a cookie has been saved to local storage. If yes then it
// triggers some action. It acts as an intermediate layer to pass communication between the datatables grid and
// the viewer.On another tab, user-selected checkboxes are saved as a cookie. Then the cookie is picked up here,
// the javacript code reads the selected values and triggers the appropriate action
var cookieMonitor = function () {
    return function () {
        var state = JSON.parse(localStorage['updated_state']);
        var enter = state.enter,
            exit = state.exit;

        exit.forEach(d => _removeOverlay(d));
        enter.forEach(d => _addOverlay(d))

    };
}();
$(window).on("storage", cookieMonitor); // via jQuery. can be used also with VanillaJS


function legendControl() {
    var legendLink = document.querySelector(`#legend-link`);

    if (!legend_added) {
        legendLink.addEventListener(`click`, () => {
            // Opens the page and stores the opened window instance
            legendWindow = window.open(`./dashboard/my_datatable.html`); // <--- THATS NEEDS TO BE EXPOSED TO THE USER. MOVE I INSIDE config.js MAYBE
        });
    }
    legend_added = true;

    $('#legend').show()
}


function closeLegend() {
    $('#legend').hide();

    if (legendWindow) {
        legendWindow.close();

        // Preventing weird behaviors
        // legendWindow = null;
    }
}

function truncateStr(strIn){
    var n = 2;
    return myUtils().fw_stripper(strIn, n);
}

var shortNames = d3.map(classColorsCodes(), d => truncateStr(d.className))
    .keys()
    .filter(d => d != "Other")
    .sort();

shortNames.forEach((d, i) => {
    // make some pixiGraphics (aka containers) objects to hold the cell polygons and name them based on their short names
    // these are just empty right now, they only have a name
    var c = new PIXI.Graphics();
    c.name = d;
    container_array.push(c)
});


run();

function run() {
    console.log('app starts')
    configSettings = config().get('default');
    var boundariesjson = configSettings.cellBoundaries;
    var celljson = configSettings.cellData;
    var q = d3.queue();
        q = q.defer(d3.json, boundariesjson);
        for (var i = 0; i < 154; i++) { // DO NOT FORGET TO REMOVE THAT (154)
            q = q.defer(d3.json, configSettings.spot_json(i));
            q = q.defer(d3.json, configSettings.cell_json(i));
        }
    q.await(onCellsLoaded(configSettings));
}


function onCellsLoaded(cfg) {
    var data_1 = [],
        // all_geneData = [],
        data_3 = [];
    return (err, ...args) => {
        args.forEach((d, i) => {
            i === 0? data_1=d:
                i % 2 === 0? data_3 = [...data_3, ...d]:
                     all_geneData = [...all_geneData, ...d]
        });
        [cellBoundaries, cellData] = postLoad([data_1, data_3]);
        dapiChart(cfg);
    }
}

function postLoad(arr) {
    //Do some basic post-processing/cleaning

    var _cellBoundaries = arr[0];
    //for some reason which I havent investigated, some cells do not have boundaries. Remove those cells
    var null_coords = _cellBoundaries.filter(d => {
        return d.coords === null
    });
    if (null_coords) {
        null_coords.forEach(d => {
            console.log('Cell_id: ' + d.cell_id + ' doesnt have boundaries')
        })
    }

    // If you need to remove the cells with no boundaries uncomment the line below:
    // _cellBoundaries = _cellBoundaries.filter(d => { return d.coords != null });

    var _cellData = arr[1];
    //stick the aggregated metrics
    var agg = aggregate(_cellData);
    _cellData.forEach((d, i) => {
        d.topClass = d.ClassName[maxIndex(d.Prob)]; // Keeps the class with the highest probability
        d.agg = agg[i];
    });

    // make sure the array is sorted by Cell_Num
    _cellData.sort(function(a,b){return a.Cell_Num-b.Cell_Num});

    return [_cellBoundaries, _cellData]
}

function maxIndex(data){
    //returns the index of the max of the input array.
    return data.reduce((iMax, x, i, arr) => x > arr[iMax] ? i : iMax, 0);
}

function aggregate(data) {
    var cellColorRamp = classColorsCodes();
    var cellColorMap = d3.map(cellColorRamp, function (d) {
        return d.className;
    });


//     for (var i = 0; i < data.length; ++i) {
//         data[i].managedData = managedData[i]
//     }
    function aggregator(data) {
        var out;
        out = d3.nest()
            .key(function (d) {
                return d.IdentifiedType;
            }) //group by IdentifiedType
            .rollup(function (leaves) {
                return {
                    Prob: d3.sum(leaves, function (d) {
                        return d.Prob;
                    }), //sum all the values with the same IdentifiedType
                    color: leaves[0].color //Get the first color code. All codes with the same IdentifiedType are the same anyway
                }
            }).entries(data)
            .map(function (d) {
                return {IdentifiedType: d.key, Prob: d.value.Prob, color: d.value.color};
            });

        // sort in decreasing order
        out.sort(function (x, y) {
            return d3.ascending(y.Prob, x.Prob);
        });

        return out
    }

    function dataManager(data) {
        var chartData = [];
        for (var i = 0; i < data.length; ++i) {
            var temp = [];
            for (var j = 0; j < data[i].ClassName.length; ++j) {
                // console.log(data[i].ClassName[j])
                temp.push({
                    IdentifiedType: cellColorMap.get(data[i].ClassName[j]).IdentifiedType,
                    color: cellColorMap.get(data[i].ClassName[j]).color,
                    Prob: data[i].Prob[j] ? data[i].Prob[j] : [data[i].Prob] //Maybe that one is better
                })
            }
            var agg = aggregator(temp);
            chartData.push({
                X: data[i].X,
                Y: data[i].Y,
                GeneCountTotal: data[i].CellGeneCount.reduce((a, b) => a + b, 0), //get the sum of all the elements in the array
                IdentifiedType: agg[0].IdentifiedType,
                color: agg[0].color,
                Prob: agg[0].Prob,
                // renderOrder: sectionFeatures.renderOrder(agg[0].IdentifiedType),

            })
        }

        return chartData
    }

    var aggData = dataManager(data);
    return aggData
}
