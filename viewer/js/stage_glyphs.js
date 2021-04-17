function _removeOverlay(name) {
    if (geneOverlays) {
        var el = geneOverlays.filter(d => d.geneName === name);
        if (el) {
            geneLayers.removeLayer(el[0].geneLayer)
        }
    }
}

function _addOverlay(name) {
    if (geneOverlays) {
        var el = geneOverlays.filter(d => d.geneName === name);
        if (el) {
            geneLayers.addLayer(el[0].geneLayer)
        }
    }
}

function refresh() {
    // if localStorage hide the relevant layers
    // Call this to ensure that when you uncheck a gene from the gene panel and then you zoom in/out
    // you will not show tha gene(s) that were not selected.
    if (localStorage['updated_state'] && JSON.parse(localStorage['updated_state']).deselected) {
        var exit = JSON.parse(localStorage['updated_state']).deselected
        exit.forEach(d => {
            _removeOverlay(d)

            if (masterMarkerContainer) {
                var x = masterMarkerContainer.getChildByName(d);
                if (x) {
                    x.visible = false
                }
            }
        })
    }
}

function renderGlyphs(evt, config) {
    var sw = evt.target.getBounds().getSouthWest();
    var ne = evt.target.getBounds().getNorthEast();
    var aw = activeWindow(sw, ne);

    console.log('filtering spots using spatial index')
    geneData = spotsIndex.range(aw.bottomLeft.x, aw.bottomLeft.y, aw.topRight.x, aw.topRight.y).map(id => all_geneData[id]);
    renderChart(geneData);
    // refresh();



    function activeWindow(sw, ne) {
        var bottomLeft = dapiConfig.t.untransform(L.point([sw.lng, sw.lat]));
        var topRight = dapiConfig.t.untransform(L.point([ne.lng, ne.lat]));

        return {'bottomLeft': bottomLeft, 'topRight': topRight}
    }



    function renderChart(geneData) {
        myDots = make_dots(geneData);
        geneData = null; //delete geneData, not needed anymore
        geneOverlays = []; // gene to layer map
        console.log('New dot layer added!')

        console.log('Adding geneLayers')
        dapiConfig.removeLayer(geneLayers);
        geneLayers = new L.LayerGroup().addTo(map);
        //create marker layer and display it on the map
        var genes = glyphSettings().map(d => d.gene).sort();
        for (var i = 0; i < genes.length; i += 1) {
            var geneName = genes[i];
            var geneLayer = L.geoJson(myDots, {
                pointToLayer: function (feature, latlng) {
                    return new svgGlyph(latlng, dapiConfig.style(feature, 'gene')).bindTooltip(feature.properties.Gene, {className: 'myCSSClass'});
                },
                filter: geneFilter(geneName),
                class: function (feature, latlng) {
                    return feature.properties.Gene + '_glyph'
                }, // i dont think i need that anymore, ok to remove (assigning a class I mean....)
                onEachFeature: onEachDot
            });
            geneLayer.addTo(geneLayers);
            geneOverlays.push({geneName, geneLayer})
        }
        console.log('geneLayers added!')
    }


    function onEachDot(feature, layer) {
        layer.on({
            mouseover: glyphMouseOver,
            mouseout: glyphMouseOut,
            click: clickGlyph,
            // popupopen: onPopupOpen,
        });
        if (feature.properties.neighbour) {
            // layer.bindPopup(spotPopup, customOptions)
        }

    }

    function glyphMouseOver(e) {
        var layer = e.target;
        dotStyleHighlight = highlightStyle(layer.feature);
        layer.setStyle(dotStyleHighlight);
        console.log('glyph mouseOver')
    }


    function glyphMouseOut(e) {
        console.log('glyph mouseout')
        resetDotHighlight(e)
        // if (map.hasLayer(glyphToNeighbourLayer)){
        //     map.removeLayer(glyphToNeighbourLayer)
        // }
    }

//create highlight style, with darker color and larger radius
    function highlightStyle(feature) {
        return {
            radius: dapiConfig.getRadius(feature.properties.size) * 2,
            fillColor: "#FFCE00",
            color: "#FFCE00",
            // color:feature.properties.glyphColor,
            weight: 3,
            opacity: 1,
            fillOpacity: 0.5
        };
    }

    function resetDotHighlight(e) {
        var layer = e.target;
        dotStyleDefault = dapiConfig.style(layer.feature);
        layer.setStyle(dotStyleDefault);
    }

    function HighlightNearbyCells(neighbours) {
        if (map.hasLayer(nearbyCellLayer)) {
            map.removeLayer(nearbyCellLayer)
        }
        nearbyCellLayer = L.geoJSON(cellPolygons, {
            filter: isCellNearby(neighbours),
            style: function (feature) {
                return {
                    fillOpacity: 0.1,
                    color: feature.properties.agg.color,
                    weight: 3,
                };
            },
            interactive: false
        });
        nearbyCellLayer.addTo(map);

        function isCellNearby(neighbours) {
            return function (feature) {
                var nearbyCells = d3.map(neighbours, d => +d.Cell_Num).keys()
                    .map(el => +el)
                    .filter(d => d != cellWatch);
                return nearbyCells.includes(feature.properties.cell_id)
            }
        }
    }

    function geneFilter(geneName) {
        return function (feature) {
            return feature.properties.Gene === geneName;
        }
    }

    function clickGlyph(e) {
        console.log('glyph clicked');
        if (glyphToNeighbourLayer) {
            map.removeLayer(glyphToNeighbourLayer)
        }

        var neighbours = e.target.feature.properties.neighbours,
            focus_cell = e.target.feature.properties.neighbour,
            focus_point = e.target.feature.geometry.coordinates,
            res = [];

        // Highlight nearby cells
        HighlightNearbyCells(neighbours);


        points_arr = neighbours.forEach(n => {
            console.log(n.Cell_Num);

            var temp_cell = cellData.filter(d => d.Cell_Num === n.Cell_Num);
            if (Array.isArray(temp_cell) && temp_cell.length) {
                res.push({
                    'X': temp_cell[0].X,
                    'Y': temp_cell[0].Y,
                    'color': temp_cell[0].agg.color,
                    'Cell_Num': temp_cell[0].Cell_Num
                })
            }
        });

        var glyphToNeighbourFc = glyphToNeighbour(res, focus_point, neighbours);
        glyphToNeighbourLayer = L.geoJson(glyphToNeighbourFc, {
            onEachFeature: onEachGlyphToNeighbourLine,
            interactive: false,
        });

        function onEachGlyphToNeighbourLine(feature, layer) {
            layer.setText('Prob = ' + d3.format("0.0%")(feature.properties.Prob), {
                offset: -5,
                center: true,
                orientation: feature.properties.angle,
                attributes: {fill: feature.properties.color, 'font-size': 16},
            });
            layer.setStyle({
                'color': feature.properties.color,
                'weight': feature.properties.Prob / 0.25 + 1,
                'className': 'geneToCellLine',
            });
        }

        glyphToNeighbourLayer.addTo(map);
        console.log('lines drawn');


        // $('.geneToCellLine')
        //     .append("text")
        //         .attr("x", (50 / 2))
        //         .attr("y", 50) //set your y attribute here
        //         .style("text-anchor", "middle")
        //         .style("font-size", "10px")
        //         .style('fill', '#707070')
        //         .text("Select an area by dragging across the lower chart from its edges");

        // stop the click from propagating any further. Mouse clicks
        // make the two controls at the bottom of the map sticky. If you want them to
        // remain frozen when you click a glyph, keep the line below
        e.stopPropagation(); //<---- that DOES NOT work, gives an error


        // i dont remember if I need these 2 lines below ( i think i dont...)
        // neighbours.forEach(d => cellData.filter(d => d.Cell_Num === focus_id))
        // // cellData.filter(d => d.Cell_Num === e.target.feature.properties.neighbours[0].Cell_Num)
    }

    function glyphToNeighbour(destination, origin, neighbours) {
        var out = {
            type: "FeatureCollection",
            features: []
        };
        for (var i = 0; i < destination.length; ++i) {
            var x = +destination[i].X,
                y = +destination[i].Y,
                lp = dapiConfig.t.transform(L.point([x, y])),
                toPoint = [lp.x, lp.y];
            var g = {
                "type": "LineString",
                "coordinates": [origin, toPoint]
            };

            //create feature properties
            var p = {
                "fromPoint": origin,
                "toPoint": toPoint,
                "color": destination[i].color,
                "Prob": d3.map(neighbours, d => d.Cell_Num).get(destination[i].Cell_Num).Prob,
                "angle": getAngle(toPoint[0], toPoint[1], origin[0], origin[1]),
            };

            console.log('origin: ' + origin);
            console.log('toPoint: ' + toPoint);
            console.log('angle: ' + p.angle);
            console.log('prob: ' + p.Prob);
            console.log("-----")
            //create features with proper geojson structure
            out.features.push({
                "geometry": g,
                "type": "Feature",
                "properties": p
            });
        }
        return out;
    }


    function getAngle(ax, ay, bx, by) {
        var dy = ay - by,
            dx = ax - bx;

        return dx < 0 ? 180 - getAngleDeg(dx, dy) : 360 - getAngleDeg(dx, dy)
    }


    /*
    * Calculates the angle between AB and the X axis
    * A and B are points (ax,ay) and (bx,by)
    */
    function getAngleDeg(dx, dy) {
        // var dy = ay-by,
        //     dx = ax-bx;
        var angleRad = Math.atan(dy / dx);
        var angleDeg = angleRad * 180 / Math.PI;

        return angleDeg
    }


    function make_dots(data) {
        var dots = {
            type: "FeatureCollection",
            features: []
        };
        for (var i = 0; i < data.length; ++i) {
            var x = data[i].x,
                y = data[i].y,
                gene = data[i].Gene;
            //console.log(gene)
            var lp = dapiConfig.t.transform(L.point([x, y]));
            var g = {
                "type": "Point",
                "coordinates": [lp.x, lp.y]
            };

            //create feature properties
            var p = {
                "id": i,
                "x": x,
                "y": y,
                "Gene": gene,
                "taxonomy": dapiConfig.getTaxonomy(gene),
                "glyphName": dapiConfig.getGlyphName(gene),
                "glyphColor": dapiConfig.getColor(gene),
                "block_id": data[i].block_id,
                "size": 30,
                "type": 'gene',
                "neighbour": parseFloat(data[i].neighbour), // why you are using parseFloat?? Cant remember why I did this!
                "neighbours": dapiConfig.getNeighbours(data[i].neighbour_array, data[i].neighbour_prob),
            };

            //create features with proper geojson structure
            dots.features.push({
                "geometry": g,
                "type": "Feature",
                "properties": p
            });
        }
        return dots;
    }


}