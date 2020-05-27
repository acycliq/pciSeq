function app(all_geneData, map) {
    loader.load(function (loader, resources) {
        var textures = [resources.plane.texture, resources.circle.texture, resources.bicycle.texture];
        var focusTextures = [resources.focusPlane.texture, resources.focusCircle.texture, resources.focusBicycle.texture];
        var markers = Array(2000000).fill().map((e, i) => ({
            "latitude": getRandom(0, 150000),
            "longitude": getRandom(0, 150000)
        }));

        // var legend = document.querySelector('div.legend.geometry');
        // var legendContent = legend.querySelector('.content');
        function mypixiLayer() {
            masterMarkerContainer = new PIXI.Graphics();
            var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;
            var _mycallBack = myCallback(markers, textures, focusTextures, all_geneData);
            return L.pixiOverlay(_mycallBack, masterMarkerContainer, {
                doubleBuffering: true,
                destroyInteractionManager: true
            });
        };

        mypixiLayer().addTo(map);
    });


    // Make an array and populate it with empty PIXI.ParticleContainers. Later these particle Containers
    // will hold the markers/spots for each gene
    if (geneContainer_array.length === 0) { // <=== TEMO WORKAROUND TO AVOID  POPULATING THIS OVER AND OVER AGAIN
        var geneNames = glyphSettings().map(d => d.gene).sort();
        geneNames.forEach((d, i) => {
            var n = all_geneData.filter(el => el.Gene === d).length;
            var c = new PIXI.ParticleContainer(n, {vertices: true});
            c.name = d;
            c.tint = markerColor(d);
            geneContainer_array.push(c)
        });
    }

    function markerColor(geneName) {
        var colorCode = glyphColor(glyphSettings().filter(d => d.gene === geneName)[0].taxonomy);
        var out = myUtils().string2hex(colorCode);
        return out
    }

    function scaleRamp(z) {
        return z === 0 ? 0.125 :
            z === 1 ? 0.125 :
                z === 2 ? 0.125 :
                    z === 3 ? 0.125 :
                        z === 4 ? 0.125 :
                            z === 5 ? 0.125 :
                                z === 6 ? 0.0625 :
                                    z === 7 ? 0.03125 : 1
    }


    function myCallback(markers, textures, focusTextures, all_geneData) {
        var firstDraw = true;
        var prevZoom;
        var markerSprites = [];

        var frame = null;
        var focus = null;
        return function (utils, event) {
            var zoom = utils.getMap().getZoom();
            var container = utils.getContainer(); // <----- I think this is referencing the same obj as '''masterMarkerContainer'''
            masterMarkerRenderer = utils.getRenderer(); // assign the renderer to the global '''masterMarkerRenderer''' variable
            var project = utils.latLngToLayerPoint;
            var scale = utils.getScale();
            var invScale = 1 / scale;
            console.log('zoom is: ' + zoom);
            console.log('scale is: ' + scale);
            console.log('inv scale is: ' + invScale);

            // get the particle containers from the geneContainer_array and add them to the container
            geneContainer_array.forEach(d => container.addChild(d));

            // get all the gene names. Now that could potentially be a problem. Maybe I should have written
            // a function just to do that. The line below gets the gene names from the particle containers name
            // which is based on the glyphConfig.js file. This file can contain a lot more (but not less!) genes that what
            // we are working with. It is a superset of our gene panel. Alternatively I could query the '''geneData''' object
            // and get my gene names but I think this is more expensive. Usually however the glyphConfig.js contains settings
            // for just our particular gene panel, however this is not necessary and is not enforced.
            var all_genes = geneContainer_array.map(d => d.name).sort();
            if (firstDraw) {
                prevZoom = zoom;

                // for each gene
                all_genes.forEach(gene => {
                    // get the relevant particle container
                    var pixiParticleContainer = geneContainer_array.filter(d => d.name === gene)[0];
                    // and filter out data not relevant to that gene
                    var tempData = all_geneData.filter(d => d.Gene === gene);
                    // Loop inside that data slice
                    tempData.forEach(function (marker) {
                        // get the coords of each datapoint and make a sprite.
                        // Then add the marker to the particle container and repeat unti exhaustion
                        // Note: The marker has to be added to the relevant particle container.
                        // DO NOT ADD THE MARKER DIRECTLY TO THE CONTAINER (aka '''masterMarkerContainer''')
                        var p = dapiConfig.t.transform(L.point([marker.x, marker.y]));
                        var coords = project([p.y, p.x]);
                        var markerSprite = new PIXI.Sprite(textures[0]);
                        markerSprite.x0 = coords.x;
                        markerSprite.y0 = coords.y;
                        markerSprite.x = coords.x;
                        markerSprite.y = coords.y;
                        markerSprite.anchor.set(0.5, 0.5);
                        pixiParticleContainer.addChild(markerSprite); //<===== HERE THE MARKER IS ADDED TO THE PARTICLECONTAINER
                        markerSprites.push(markerSprite);
                        markerSprite.legend = marker.city || marker.label;
                    });
                })
                var quadTrees = {};
                var start = +new Date();
                for (var z = map.getMinZoom(); z <= map.getMaxZoom(); z++) {
                    if (z < 9){
                        quadTrees[z] = null  // I do not need the tree, I wont do mouse events when zoom < 7
                    }
                    else {
                        var rInit = ((z <= 7) ? 10 : 24) / utils.getScale(z);
                        quadTrees[z] = solveCollision(markerSprites, {r0: rInit, zoom: z});
                    }
                }
                var end = +new Date();
                var time = end - start;
                console.log('total execution time = ' + time / 1000 + 'sec');

                function findMarker(ll) {
                    var layerPoint = project(ll);
                    var quadTree = quadTrees[utils.getMap().getZoom()];
                    var marker;
                    var rMax = quadTree.rMax;
                    var found = false;
                    quadTree.visit(function (quad, x1, y1, x2, y2) {
                        if (!quad.length) {
                            var dx = quad.data.x - layerPoint.x;
                            var dy = quad.data.y - layerPoint.y;
                            var r = quad.data.scale.x * 16;
                            if (dx * dx + dy * dy <= r * r) {
                                marker = quad.data;
                                found = true;
                            }
                        }
                        return found || x1 > layerPoint.x + rMax || x2 + rMax < layerPoint.x || y1 > layerPoint.y + rMax || y2 + rMax < layerPoint.y;
                    });
                    return marker;
                }

                map.on('click', function (e) {
                    var redraw = false;
                    if (focus) {
                        focus.texture = textures[focus.textureIndex];
                        focus = null;
                        L.DomUtil.addClass(legend, 'hide');
                        legendContent.innerHTML = '';
                        redraw = true;
                    }
                    var marker = findMarker(e.latlng);
                    if (marker) {
                        marker.texture = focusTextures[marker.textureIndex];
                        focus = marker;
                        legendContent.innerHTML = marker.legend;
                        L.DomUtil.removeClass(legend, 'hide');
                        redraw = true;
                    }
                    if (redraw) utils.getRenderer().render(container);
                });

                var self = this;
                map.on('mousemove', L.Util.throttle(function (e) {
                    var marker = findMarker(e.latlng);
                    if (marker) {
                        L.DomUtil.addClass(self._container, 'leaflet-interactive');
                    } else {
                        L.DomUtil.removeClass(self._container, 'leaflet-interactive');
                    }
                }, 32));
            }

            var start = null;
            var delta = 250;

            function animate(timestamp) {
                var progress;
                if (start === null) start = timestamp;
                progress = timestamp - start;
                var lambda = progress / delta;
                if (lambda > 1) lambda = 1;
                lambda = lambda * (0.4 + lambda * (2.2 + lambda * -1.6));
                markerSprites.forEach(function (markerSprite) {
                    var targetScale = zoom <= 7 ? scaleRamp(zoom) : 1 / (2 * utils.getScale(event.zoom));
                    markerSprite.scale.set(targetScale);
                });
                masterMarkerRenderer.render(container);
                if (progress < delta) {
                    frame = requestAnimationFrame(animate);
                }
            }

            if (firstDraw || prevZoom !== zoom) {
                markerSprites.forEach(function (markerSprite) {
                    var targetScale = zoom <= 7 ? scaleRamp(zoom) : 1 / (2 * utils.getScale(event.zoom));
                    markerSprite.scale.set(targetScale);
                });
            }

            // Not quite sure if that makes any difference, maybe I am not using it right
            // if (!firstDraw && prevZoom !== zoom) {
            //     frame = requestAnimationFrame(animate);
            // }


            firstDraw = false;
            prevZoom = zoom;
            masterMarkerRenderer.render(container);
        }
    }
};


// myMarkersContainer = masterMarkerContainer.getChildByName('myMarkers');
// myMarkersContainer.visible = false
// masterMarkerRenderer.render(masterMarkerContainer)


// [...new Set(all_geneData.map(d => d.Gene))].sort()