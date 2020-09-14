function drawPoly(color, alpha, project, mycontainer) {
    if(isNaN(color)){
        console.log('Oh!')
    }
    // mycontainer is now an empty pixiraphics object, it has only a name, the shortnened cell type name.
    return function (poly) {
        var shape = new PIXI.Polygon(poly[0].map(function (point) {
            var proj = project([point[1], point[0]]);
            return new PIXI.Point(proj.x, proj.y);
        }));
        mycontainer.beginFill(color, alpha);
        mycontainer.drawShape(shape);
        if (poly.length > 1) {
            for (var i = 1; i < poly.length; i++) {
                var hole = new PIXI.Polygon(poly[i].map(function (point) {
                    var proj = project([point[1], point[0]]);
                    return new PIXI.Point(proj.x, proj.y);
                }));
                mycontainer.drawShape(hole);
                mycontainer.addHole();
            }
        }
    };
}

function drawCellOutline(color, alpha, project, container, coords) {
    // Draws the outlines/boundaries of the cell, but I dont think I am going to be using it, there is no need
    container.lineStyle(2, color, 1);
    coords[0].forEach(function (point, index) {
        var proj = project([point[1], point[0]]);
        if (index === 0) container.moveTo(proj.x, proj.y);
        else container.lineTo(proj.x, proj.y);
    });
}

function drawCellPolygons() {
    cellPolygons = myUtils().poly_collection(cellBoundaries, dapiConfig.t);
    cellBoundaries = null;
    masterCellContainer = new PIXI.Graphics(); // Assign this to the global variable 'masterCellContainer'
    var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;

    var poly = renderPolygons(cellPolygons);
    polygonsOverlay = L.pixiOverlay(poly, masterCellContainer, {
        doubleBuffering: doubleBuffering,
        destroyInteractionManager: true
    });

    // polygonsOverlay.addTo(map);
    return polygonsOverlay
}


function renderPolygons(data) {
    console.log('In renderPolygons')
    var alphaScale = d3.scaleLinear()
        .domain([0, map.getMinZoom(), map.getMaxZoom()])
        .range([1, 1, 0.2]); // <-- huh?? I dont remember why i did this
    alphaScale.clamp(true);
    var firstDraw = true;

    return function (utils) {
        console.log('you passed in ' + data.features.length + ' data points')
        var container = utils.getContainer();  // That has to be pointing to the same object as the global var 'pixiContainer' (hmm, masterCellContainer I probably mean! Not pixiContainer)
        masterCellRenderer = utils.getRenderer();   // Assign this to the global variable 'masterCellRenderer'
        var gl = masterCellRenderer.gl;
        var project = utils.latLngToLayerPoint;
        var zoom = utils.getMap().getZoom();
        var scale = utils.getScale;
        var bounds;
        var tree = new RBush();
        var mouseTarget;

        // ok, thats quite crucial. I think only one pixiRenderer should be used. Every other PixiGraphics object should
        // be attached to pixiRenderer and then you can add pixiRenderer as an overlay to the leaflet map
        cellContainer_array.forEach(d => {
            container.addChild(d)
            // dapiConfig.customControl._addButton(d.name) // update the layer control on the map
        });
        // dapiConfig.customControl._isEnabled = false; //prevents the control from adding new elements

        if (firstDraw) {
            if ( (masterCellRenderer.type === PIXI.RENDERER_TYPE.WEBGL) && (mapboxgl.supported({failIfMajorPerformanceCaveat: true})) ){
                document.getElementById('webgl').innerHTML = '<span class="greenDot"> </span><span> GPU acceleration: Enabled</span>'
                console.log('GPU acceleration: Enabled')
                gl.blendFunc(gl.ONE, gl.ZERO);
                document.querySelector('#webgl').style.display = 'block';
            } else {
                document.getElementById('webgl').innerHTML = '<span class="redDot blinking"> </span><span> Browser is not using the GPU</span>';
                console.log('GPU acceleration: Disabled')
                document.querySelector('#webgl').style.display = 'block';
                // document.body.removeChild(document.querySelector('#webgl'));
            }

            // var markerCoords = project(markerLatLng)
            data.features.forEach(function (feature, index) {
                var color = myUtils().string2hex(feature.properties.agg.color),
                    alpha = 0.8;
                if (feature.geometry === null) return;
                bounds = L.bounds(feature.geometry.coordinates[0]);
                if (feature.geometry.type === 'Polygon') {
                    // 1. find the shortname of the most probable cell-type for this particular cell
                    var cName = feature.properties.topClass;

                    // 2. grab the relevant pixiGraphics object
                    var cont = cellContainer_array.filter(d => d.name === cName)[0];

                    // 3, now that you have the correct pixiGraphics object, draw the polygon.
                    // In this manner each pixiGraphics object will hold polygons that have the same shortname
                    if (cName !== "Zero"){
                        drawPoly(color, alpha, project, cont)(feature.geometry.coordinates);
                    }
                    // drawCellOutline(color, alpha, project, cont, feature.geometry.coordinates);

                } else if (feature.geometry.type === 'MultiPolygon') {
                    feature.geometry.coordinates.forEach(drawPoly(color, alpha));
                }
                tree.insert({
                    minX: bounds.min.x,
                    minY: bounds.min.y,
                    maxX: bounds.max.x,
                    maxY: bounds.max.y,
                    feature: feature,
                })
            });

            var focus,
                spots,
                lineStrings,
                lastVisited,
                clickedCell ;

            utils.getMap().on('click', onclick);

            function onclick(e){
                // pinnedControls =! pinnedControls;
                togglePin(mouseTarget);
                if (mouseTarget){
                    console.log('mouse click in some cell');
                    if(clickedCell) {map.removeLayer(clickedCell)}
                    clickedCell = highlightOutline(mouseTarget, ['#ffe513', 5])
                    clickedCell.addTo(map)
                }
                console.log('Controls are pinned: ' + pinnedControls)
            }

            utils.getMap().on('mousemove', L.Util.throttle(function (e) {
                onMousemove(e);
            }, 32));

            function togglePin(mouseTarget){
                // you can unpin from anywhere on the map
                if (pinnedControls){
                     pinnedControls = false;
                     if(clickedCell) {map.removeLayer(clickedCell)}
                     d3.selectAll('.ribbon').nodes().map(d => d.src = 'https://cdn2.iconfinder.com/data/icons/snipicons/500/pin-128.png')
                }
                else{
                    if (mouseTarget){
                        pinnedControls = true;
                        d3.selectAll('.ribbon').nodes().map(d => d.src = 'https://cdn2.iconfinder.com/data/icons/oxygen/48x48/actions/note2.png')
                    }
                }
            }

            function onMousemove(e) {
                // 1. show the mouse coords control
                $('.uiElement.label').show();

                mouseTarget = findFeature(e.latlng);
                    // isChecked =
                if (mouseTarget && myTreeControl._getSelected().includes(mouseTarget.properties.topClass)) {
                    lastVisited = mouseTarget.properties.cell_id; // <-- lastVisited the same as cellWatch, right? Why not just use cellWatch then?

                    // 1. First draw/highlight the cell boundaries
                    if (focus) map.removeLayer(focus);
                    focus = highlightOutline(mouseTarget);
                    focus.addTo(map);

                    if (!pinnedControls && !hiddenControls){
                        // show the info control, then update
                        $('.leaflet-bottom.leaflet-left').show();
                        dapiConfig.info.update(mouseTarget.properties);

                        // show the donut control, then update
                        $('.leaflet-bottom.leaflet-right').show();
                        donutchart(mouseTarget.properties);
                        renderDataTable(mouseTarget.properties);
                         // d3.select('.ribbon').node().src = 'https://cdn2.iconfinder.com/data/icons/oxygen/48x48/actions/note2.png'

                        // finally show the checkbox that hides/unhides them
                        $('.panelsToggle').show();
                    }
                    else {
                        // d3.select('.ribbon').node().src = 'https://cdn2.iconfinder.com/data/icons/snipicons/500/pin-128.png'
                    }


                    if (myDots) {
                        // 2. If you the viewer has also plotted the spots, then plot again those
                        //    few spots that are linked to the target cell
                        if (spots) map.removeLayer(spots);
                        spots = drawTargetSpots(myDots, mouseTarget);
                        spots.addTo(map);

                        // 3. Draw now the line that connect the spots to the parent cell
                        var centroid = mouseTarget.properties.centroid,
                            focusSpots;

                        // This is inefficient. You apply almost the same filter to get the spots
                        focusSpots = myDots.features.filter(d => d.properties.neighbour === mouseTarget.properties.cell_id);
                        focusSpots = dapiConfig.makeLineStringFeatures(focusSpots, centroid);
                        if (lineStrings) map.removeLayer(lineStrings);
                        lineStrings = drawLines(focusSpots);
                        lineStrings.addTo(map)

                    }
                } else {
                    if (focus) map.removeLayer(focus);
                    if (spots) map.removeLayer(spots);
                    if (lineStrings) map.removeLayer(lineStrings);
                    if (glyphToNeighbourLayer) map.removeLayer(glyphToNeighbourLayer);
                    if (nearbyCellLayer) map.removeLayer(nearbyCellLayer);

                    // reset the info control
                    if (!pinnedControls) {
                        dapiConfig.info.update();
                    }

                    // nullify the cellWatch variable
                    cellWatch = null
                }
                // console.log(feat);
                // console.log('Mouse click fired')
            }

            function highlightOutline(feature, args) {
                cellWatch = feature.properties.id; // just keep a note of the cell id which is about to be draxn
                return L.geoJSON(feature, {
                    style: function (feature) {
                        return {
                            fillOpacity: 0.1,
                            color: args? args[0]: feature.properties.agg.color,
                            weight: args? args[1]: 3,
                        };
                    },
                    interactive: false
                });
            }

            function drawTargetSpots(myDots, feat) {
                return L.geoJSON(myDots, {
                    pointToLayer: function (feature, latlng) {
                        return new svgGlyph(latlng, dapiConfig.style(feature))
                            .bindTooltip(feature.properties.Gene, {className: 'myCSSClass'})
                    },
                    filter: spotFilter(feat),
                    interactive: true
                });
            }


            function drawLines(focusSpots) {
                function onEachLineFeature(feature, layer) {
                    if (layer instanceof L.Polyline) {
                        layer.setStyle({
                            'color': feature.properties.color,
                            'weight': 1.5,
                        });
                    }
                }

                return L.geoJSON(focusSpots, {
                    onEachFeature: onEachLineFeature
                })

            }


            function spotFilter(feat) {
                // if you are still in the same cell, do not draw the glyphs again
                console.log('Filter called')
                return function (feature) {
                    if ((feature.properties.neighbour === feat.properties.cell_id)
                        & (lastVisited != feat.properties.cell_id)) {
                        return true
                    }
                }
            }

            function findFeature(latlng) {
                // var point = project(latlng);
                var point = {'x': latlng.lng, 'y': latlng.lat};
                var features = tree.search({
                    minX: point.x,
                    minY: point.y,
                    maxX: point.x,
                    maxY: point.y
                });
                for (var i = 0; i < features.length; i++) {
                    var feat = features[i].feature;
                    if (feat.geometry.type === 'Polygon') {
                        if (containsPoint(feat.geometry.coordinates, point)) return feat;
                    } else {
                        for (var j = 0; j < feat.geometry.coordinates.length; j++) {
                            var ring = feat.geometry.coordinates[j];
                            if (containsPoint(ring, point)) return feat;
                        }
                    }
                }
            }

            function containsPoint(polygon, p) {
                var inside = false,
                    part, p1, p2, i, j, k, len, len2;
                // ray casting algorithm for detecting if point is in polygon
                for (i = 0, len = polygon.length; i < len; i++) {
                    part = polygon[i];

                    for (j = 0, len2 = part.length, k = len2 - 1; j < len2; k = j++) {
                        p1 = part[j];
                        p2 = part[k];

                        if (((p1[1] > p.y) !== (p2[1] > p.y)) && (p.x < (p2[0] - p1[0]) * (p.y - p1[1]) / (p2[1] - p1[1]) + p1[0])) {
                            inside = !inside;
                        }
                    }
                }
                return inside;
            }

        }
        firstDraw = false;
        // change the opacity depending on the zoom level
        cellContainer_array.filter(d => d.alpha = alphaScale(zoom));
        // dapiConfig.customControl._refresh() // That renders the chart again, but It has to be put somewhere more prominent!
        myTreeControl._refresh()
        // pixiRenderer.render(pixiContainer);
    }
}