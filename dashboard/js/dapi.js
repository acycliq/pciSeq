function dapi(cfg) {
    console.log('Doing Dapi plot');

    var img = cfg.imageSize,
        tiles = cfg.tiles,
        roi = cfg.roi;
    // var img = [227951, 262144],
    //     roi = {"x0": 0, "x1": 40000, "y0": 0, "y1": 46000};

    var a = img[0] / (roi.x1 - roi.x0),
        b = -img[0] / (roi.x1 - roi.x0) * roi.x0,
        c = img[1] / (roi.y1 - roi.y0),
        d = -img[1] / (roi.y1 - roi.y0) * roi.y0;

    // This transformation maps a point from the roi domain to the domain defined by [0,0] amd [img[0], img[1]].
    var t = new L.Transformation(a, b, c, d);

    // The transformation in this CRS maps the the top left corner to (0,0) and the bottom right to (256, 256)
    L.CRS.MySimple = L.extend({}, L.CRS.Simple, {
        transformation: new L.Transformation(1 / 1024, 0, 1 / 1024, 0),
    });


    map = L.map('mymap', {
        crs: L.CRS.MySimple,
        attributionControl: false,
    }).setView([img[1], img[0] / 2], 2);
    L.tileLayer(tiles, {
        minZoom: 0,
        maxZoom: 10
    }).addTo(map);

    function getTaxonomy(gene) {
        if (glyphMap.get(gene)) {
            out = glyphMap.get(gene).taxonomy
        } else {
            out = glyphMap.get('Generic').taxonomy
        }
        return out
    }

    function getGlyphName(gene) {
        if (glyphMap.get(gene)) {
            out = glyphMap.get(gene).glyphName
        } else {
            out = glyphMap.get('Generic').glyphName
        }
        return out
    }

    // get the svg markers (glyphs)
    var glyphs = glyphSettings();
    var getColor = glyphColor;
    var glyphMap = d3.map(glyphs, function (d) {
        return d.gene;
    });

    //get the class colors
    var classColors = classColorsCodes();
    var classColorsMap = d3.map(classColors, function (d) {
        return d.className;
    });

    //calculate radius so that resulting circles will be proportional by area
    function getRadius(y) {
        r = Math.sqrt(y / Math.PI)
        return r;
    }

    var myRenderer = L.canvas({
        padding: 0.5,
    });

    //create style, with fillColor picked from color ramp
    function style(feature) {
        return {
            radius: getRadius(feature.properties.size),
            shape: feature.properties.glyphName,
            //fillColor: "none",//getColor(feature.properties.taxonomy),
            color: glyphColor(feature.properties.taxonomy),
            weight: 0.85,
            opacity: 0.85,
            fillOpacity: 0.0,
            renderer: myRenderer,
        };
    }

    function getNeighbours(neighbour_array, neighbour_prob) {
        var data = [];
        if (neighbour_array) {
            for (var i = 0; i < neighbour_array.length; i++) {
                data.push({
                    Cell_Num: +neighbour_array[i],
                    Prob: +neighbour_prob[i],
                })
            }

            // Sort now in decreasing order.
            data.sort(function (x, y) {
                return d3.ascending(y.Prob, x.Prob)
            })
        } else {
            data.push({
                    Cell_Num: null, // null or NaN ?
                    Prob: null,
                }
            )
        }

        return data
    }

    function removeLayer(layer) {
        if (map.hasLayer(layer)) {
            map.removeLayer(layer);
            // layer = undefined;
            console.log('Layer removed')
        }
    }

    function addLayer(layer) {
        if (map.hasLayer(layer)) {
            console.log('Already added')
        } else {
            layer.addTo(map);
            console.log('Layer added')
        }
    }

    function makeLineStringFeatures(destination, origin) {
        var o = t.transform(L.point(origin)),
            fromPoint = [o.x, o.y];
        var out = {
            type: "FeatureCollection",
            features: []
        };
        for (var i = 0; i < destination.length; ++i) {
            var spot = destination[i].properties;
            var x = +spot.x,
                y = +spot.y,
                gene = spot.Gene,
                lp = t.transform(L.point([x, y])),
                toPoint = [lp.x, lp.y];
            var g = {
                "type": "LineString",
                "coordinates": [fromPoint, toPoint]
            };

            //create feature properties
            var p = {
                "gene": gene,
                // "Cell_Num": origin.Cell_Num,
                "fromPoint": fromPoint,
                "toPoint": toPoint,
                "color": getColor(getTaxonomy(gene)),
                // "color": getColor(glyphMap.get(gene).taxonomy),
            };

            //create features with proper geojson structure
            out.features.push({
                "geometry": g,
                "type": "Feature",
                "properties": p
            });
        }
        return out;
    }


    // control that shows state info on hover
    var info = L.control({
        position: 'bottomleft'
    });

    info.onAdd = function (map) {
        this._div = L.DomUtil.create('div', 'info'); // create a div with a class "info"
        this.update();
        return this._div;
    };

    function infoMsg(cellFeatures) {
        var str1 = '</div></td></tr><tr class><td nowrap><div><b>';
        var str2 = '&nbsp </b></div></td><td><div>';
        var out = '';
        var temp = [];
        var sdata = [];
        if (cellFeatures) {
            for (var i = 0; i < cellFeatures.ClassName.length; ++i) {
                temp.push({
                    ClassName: cellFeatures.ClassName[i],
                    IdentifiedType: classColorsMap.get(cellFeatures.ClassName[i]).IdentifiedType,
                    Prob: cellFeatures.Prob[i],
                })
            }

            temp = d3.nest()
                .key(function (d) {
                    return d.IdentifiedType;
                })
                .rollup(function (leaves) {
                    return d3.sum(leaves, function (d) {
                        return d.Prob;
                    })
                }).entries(temp)
                .map(function (d) {
                    return {IdentifiedType: d.key, Prob: d.value};
                });

            // sort in decreasing order
            temp.sort(function (x, y) {
                return d3.ascending(y.Prob, x.Prob);
            })

            for (var i = 0; i < temp.length; i++) {
                out += str1 + temp[i].IdentifiedType + str2 +
                    Math.floor(temp[i].Prob * 10000) / 100 + '%'
            }
        } else {
            // do nothing
        }
        return out;
    };

    // method that we will use to update the control based on feature properties passed
    info.update = function (cellFeatures) {
        var msg = infoMsg(cellFeatures);
        this._div.innerHTML = '<div class="infoTitle"><h4>Cell Info</h4><img src="https://cdn2.iconfinder.com/data/icons/snipicons/500/pin-128.png" class="ribbon"></div>' + (cellFeatures ?
                '<table style="width:110px;">' +
                '<tbody><tr style="width:110px; border-bottom:1px solid Black; font-weight: bold"><td><div><b>Class </b></div></td><td><div> Prob' +
                msg +
                '<tbody><tr style="width:110px; border-top:1px solid black;"><td><div><b>Cell Num: </b></div></td><td><div>' + cellFeatures.cell_id +
                '</div></td></tr></tbody></table>' :
                '<b>Hover over  a cell</b>'

        );
    };

    // control that shows donut chart on hover
    var summary = L.control({
        position: 'bottomright'
    });

    summary.onAdd = function (map) {
        var divhtml = dapiConfig.createDiv()
        var container = document.createElement('div');
        container.setAttribute("id", "container");
        container.innerHTML = divhtml;
        this._div = container.firstChild;
        this.update();
        return this._div;
    };


    // control that shows datatable chart on hover
    var datatable = L.control({
        position: 'bottomright'
    });

    datatable.onAdd = function (map) {
        this._div = L.DomUtil.create('table', 'display compact custom');
        this._div.id = "dtTable"
        d3.select(this._div).attr('data-page-length', '5')
        // d3.select(this._div).append('svg').attr('width', '280').attr('height', '200')
        this.update();
        return this._div;
    };

    // method that we will use to update the control based on feature properties passed
    summary.update = function (x) {
        if (x) {
            this._div.innerHTML = x.innerHTML
        }
    };

    // method that we will use to update the control based on feature properties passed
    datatable.update = function (x) {
        if (x) {
            this._div.innerHTML = x.innerHTML
        }
    };


    function createDiv() {
        // This will hold the two charts on the bottom right corner
        // lots of styling in here. Better keep these in the main html of css file
        var myDiv = "<div class='tab-pane active fade in' id='map-content'>" +
                        "<div class='container-fluid'>" +
                            "<div class='col-sm-12'>" +
                                "<div class='row'>" +
                                    "<div class='myTitle' id='dtTitle' style='margin-bottom:5px'> <h4>Highlighted Cell</h4>  " +
                                    " <img src='https://cdn2.iconfinder.com/data/icons/snipicons/500/pin-128.png' class='ribbon'/> " +
                                "</div>" +
                                "</div>" +
                                "<div class='row' style='background-color: rgba(255, 255, 255, 0.0)'>" +
                                    "<div class='chart-wrapper'>" +
                                    // "<div class='chart-title' id='dtTitle'> </div>" +
                                        "<div class='chart-stage'>" +
                                            "<div class='col-sm-5'>" +
                                                "<div class='chart-stage' style='background-color: rgba(255, 255, 255, 0.8)'>" +
                                                    "<table id='dtTable' class='display compact custom' data-page-length='5' width=100%></table>" +
                                                "</div>" +
                                            "</div>" +
                                            "<div class='col-sm-7'>" +
                                                "<div class='chart-stage' style='background-color: rgba(255, 255, 255, 0.0)'> " +
                                                    "<div class='summary' id='pie'> " +
                                                        "<svg width='300' height='180'></svg>" +
                                                    "</div>" +
                                                "</div>" +
                                            "</div>" +
                                        "</div>" +
                                    "</div>" +
                                "</div>" +
                            "</div>" +
                        "</div>" +
                    "</div>";

        return myDiv
    }



    var toggleMapControl = L.control({
        position: 'topright'
    });


    toggleMapControl.onAdd = function (map) {
        var div = L.DomUtil.create('div');
        div.innerHTML =
        '<div class="leaflet-control-layers leaflet-control-layers-expanded"> ' +
        '  <form> ' +
        '    <input class="leaflet-control-layers-overlays" id="command"  ' +
        '      onclick = dapiConfig.toggleMapControl.update(this.checked) type="checkbox"> ' +
        '      Hide Dapi ' +
        '    </input> ' +
        '  </form> ' +
        ' </div>';
        return div;
    };

    toggleMapControl.update = function (bool) {
        if (bool) {
            $('.leaflet-tile-container').hide();
            console.log('Background image: hidden')
        } else {
            $('.leaflet-tile-container').show();
            console.log('Background image: visible')
        }

    };


    function tree(data) {
        // makes the tree object to pass into the tree control as an overlay
        var mapper = {},
            root = {
                label: 'Cell Classes',
                selectAllCheckbox: 'Un/select all',
                children: []
            };

        for (var str of data) {
            var sep = configSettings.class_name_separator,
                splits,
                label = '';
            // let splits = str.match(/[a-zA-Z]+|[0-9]+/g), //str.split('.'),
            if (sep === '') {
                console.log('Assuming that class name is a string followed by a number, like Astro1, Astro2 etc');
                splits = str.match(/[a-zA-Z]+|[0-9]+/g) //str.split('.'),
            }
            else {
                splits = str.split(sep)
            };
            splits.reduce(myReducer(label), root)
        }

        function myReducer(label) {
            return function (parent, place, i, arr) {
                if (label) {
                    var sep = configSettings.class_name_separator;
                    label += sep + `${place}`; // `.${place}`;
                }
                else
                    label = place;

                if (!mapper[label]) {
                    var o = {label: label};
                    o.collapsed = true;
                    if (i === arr.length - 1) {
                        o.layer = cellContainer_array.filter(d=> d.name === label)[0]
                        // o.layer = masterCellContainer.getChildByName(label);
                    }
                    mapper[label] = o;
                    parent.selectAllCheckbox = true;
                    parent.children = parent.children || [];
                    parent.children.push(o)
                }
                return mapper[label];
            }
        }

        return root
    }


     function treeControl(data) {
        return L.control.layers.tree({}, tree(data), {position:'topright'});
     }



    // add the customised control
    // customControl = L.control.custom().addTo(map);

    var dapiData = {};
    dapiData.map = map;
    dapiData.t = t;
    dapiData.style = style;
    dapiData.getTaxonomy = getTaxonomy;
    dapiData.getGlyphName = getGlyphName;
    dapiData.getNeighbours = getNeighbours;
    dapiData.getColor = getColor;
    dapiData.getRadius = getRadius;
    dapiData.removeLayer = removeLayer;
    dapiData.addLayer = addLayer;
    dapiData.makeLineStringFeatures = makeLineStringFeatures;
    dapiData.info = info;
    dapiData.summary = summary;
    dapiData.toggleMapControl = toggleMapControl;
    dapiData.createDiv = createDiv;
    dapiData.datatable = datatable;
    // dapiData.customControl = customControl;
    dapiData.treeControl = treeControl;
    return dapiData
}

function dapiChart(config) {

    dapiConfig = dapi(config);
    var map = dapiConfig.map;
    map.whenReady(d => {
        console.log('Map is ready')
    });
    map.on('moveend', moveend(config)); // Show the polygons when move has ended
    map.on('movestart', movestart); // Hide the polygons when you are about to move the chart. Thats a trick to make the chart look lighter when the number pf polygons is toooo large
    map.on('zoomanim', zoomanim_end);
    map.on('zoomend', zoomanim_end); // attach zooanim_end to both zoomanim and zoomend

    // Now add the info control  to map...
    dapiConfig.info.addTo(map);

    //...and the summary control too
    dapiConfig.summary.addTo(map);

    // and the toggle to hide/show the background image
    // dapiConfig.toggleMapControl.addTo(map); // changed my mind. I dont like the way it is placed, I did another button for this, simple one, not L.control

    //... and show the legend button
    legendControl();

    //... and show the button to hide/show the dapi and the pie/info panels
    $('#hideDapiAndPanels').show();
    console.log('check boxes added');

    cellClasses = [...new Set(cellData.map(d => d.topClass))].sort();
    cellClasses.forEach((d, i) => {
        // make some pixiGraphics (aka containers) objects to hold the cell polygons and name them based on their short names
        // these are just empty right now, they only have a name
        var c = new PIXI.Graphics();
        c.name = d;
        c.options = [];
        c.options.minZoom = 0;  // Needed only to fool the layer tree control and prevent an error from being raised
        c.options.maxZoom = 10; // Needed only to fool the layer tree control and prevent an error from being raised
        cellContainer_array.push(c)
    });

    // that needs to be created before we do the cell polygons (Should do that in a better way, maybe inside drawCellPolugons?)
    // Add the control to switch on/off the cell polygons per class
    myTreeControl = dapiConfig.treeControl(cellClasses);
    myTreeControl.addTo(map);
    myTreeControl._checkAll();

    // 1. draw the cell polygons
    cellPolyLayer = drawCellPolygons();
    cellPolyLayer.addTo(map);
    console.log('cellPolyLayer added to the map');


    // draw the spots
    // add_spots(all_geneData, map);

    // draw the spots (particle Containers approach)
    add_spots_patched(all_geneData, map);


    function moveend(config) {
        console.log('Triggering moveend callback');

        if (masterCellContainer){
            console.log('setting cell container to true')
            masterCellContainer.children.map(d => d.visible = true)
        }

        // if zoom >= 7 then render the glyphs too
        return function (evt) {
            if (map.getZoom() >= zoomSwitch) {
                // hide the markers drawn by pixi
                geneContainer_array.map(d => d.visible = false);

                //and then render the glyphs (leaflet + canvas)
                renderGlyphs(evt, config);
                refresh();

            } else {
                dapiConfig.removeLayer(geneLayers);
                // closeLegend()
                // localStorage.clear();

                // show the markers drawn by pixi
                geneContainer_array.map(d => d.visible = true);

                // call refresh(). If you have unchecked a gene from the gene panel
                // then the spots from that gene should not show up
                refresh()
            }
            console.log("Current Zoom Level =" + map.getZoom());
            console.log('exiting moveend callback');
            console.log('')
        };
    }

    function movestart(evt){
        console.log('map was moved.')
        console.log('setting cell container to false');
        masterCellContainer.children.map(d => d.visible = false)
    }

    function zoomanim_end(){
        // make sure dapi remains hidden when you change zoom levels and the 'Hide Dapi' checkbox is checked
        dapiConfig.toggleMapControl.update( document.getElementById('dapiToggle').checked );
    }


    // make placeholder for the coordinates control
    addControlPlaceholders(map);

    function addControlPlaceholders(map) {
        var corners = map._controlCorners,
            l = 'leaflet-',
            container = map._controlContainer;

        function createCorner(vSide, hSide) {
            var className = l + vSide + ' ' + l + hSide;
            corners[vSide + hSide] = L.DomUtil.create('div', className, container);
        }

        createCorner('verticalcenter', 'left');
        createCorner('verticalcenter', 'right');
        createCorner('verticalcenter', 'horizontalcenter');
        createCorner('bottom', 'horizontalcenter');
        createCorner('top', 'horizontalcenter');

    }

    // Patch first to avoid longitude wrapping.
    L.Control.Coordinates.include({
        _update: function (evt) {
            var pos = evt.latlng,
                opts = this.options;
            if (pos) {
                //pos = pos.wrap(); // Remove that instruction.
                // Get the mouse location in actual coordinates and then express it as a latLng object
                var coords = dapiConfig.t.untransform(L.point([pos.lng, pos.lat])),
                    pos = L.latLng(coords.y, coords.x)
                this._currentPos = pos;
                this._inputY.value = L.NumberFormatter.round(pos.lng, opts.decimals, opts.decimalSeperator);
                this._inputX.value = L.NumberFormatter.round(pos.lat, opts.decimals, opts.decimalSeperator);
                this._label.innerHTML = this._createCoordinateLabel(pos);
            }
        }
    });

    L.control.coordinates({
        position: "tophorizontalcenter", //optional default "bootomright"
        // position: "topright", //optional default "bootomright"
        decimals: 0, //optional default 4
        decimalSeperator: ".", //optional default "."
        labelTemplateLat: "y: {y}", //optional default "Lat: {y}"
        labelTemplateLng: "x: {x}", //optional default "Lng: {x}"
        enableUserInput: true, //optional default true
        useDMS: false, //optional default false
        useLatLngOrder: false, //ordering of labels, default false-> lng-lat
        markerType: L.marker, //optional default L.marker
        markerProps: {}, //optional default {},
        // labelFormatterLng : function(lng){return lng+" lng"}, //optional default none,
        // labelFormatterLat : function(lat){return lat+" lat"}, //optional default none
        // customLabelFcn: function(latLonObj, opts) { "Geohash: " + encodeGeoHash(latLonObj.lat, latLonObj.lng)} //optional default none
    }).addTo(map);

    // Hide the controls the first time the page loads up
    $('.uiElement.label').hide();
    // $('#legend').hide();
    $('.leaflet-bottom.leaflet-left').hide();
    $('.leaflet-bottom.leaflet-right').hide();
    $('.panelsToggle').hide()


}