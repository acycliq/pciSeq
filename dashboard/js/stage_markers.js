function add_spots(all_geneData, map) {

    loader.load(function (loader, resources) {
        var textures = [resources.plane.texture, resources.circle.texture, resources.bicycle.texture];
        var focusTextures = [resources.focusPlane.texture, resources.focusCircle.texture, resources.focusBicycle.texture];

        // var legend = document.querySelector('div.legend.geometry');
        // var legendContent = legend.querySelector('.content');
        function mypixiLayer() {
            masterMarkerContainer = new PIXI.Graphics();
            var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;
            var _mycallBack = myCallback(textures, focusTextures, all_geneData);
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
            var c = new PIXI.particles.ParticleContainer(n, {vertices: true});
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


    function myCallback(textures, focusTextures, all_geneData) {
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
            var tree = new RBush(); // <--- KDBush looks like a better option maybe!
            var scale = utils.getScale();
            var invScale = 1 / scale;
            console.log('zoom is: ' + zoom);
            console.log('scale is: ' + scale);
            console.log('inv scale is: ' + invScale);


            // var p = new PIXI.Graphics();
            // p.beginFill(0x000000);
            // p.lineStyle(0);
            // p.drawCircle(100, 100, 10);
            // p.endFill();
            //
            // // var t = PIXI.RenderTexture.create(p.width, p.height);
            // // masterMarkerRenderer.render(p, t);
            // var t = masterMarkerRenderer.generateTexture(p);
            // var myNewSprite = new PIXI.Sprite(t);


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
                        var markerSprite = new PIXI.Sprite(textures[1]);
                        markerSprite.x = coords.x;
                        markerSprite.y = coords.y;
                        markerSprite.anchor.set(0.5, 0.5);
                        pixiParticleContainer.addChild(markerSprite); //<===== HERE THE MARKER IS ADDED TO THE PARTICLECONTAINER
                        markerSprites.push(markerSprite);
                        markerSprite.legend = marker.city || marker.label;

                        // // comment this in for now. Will get back later
                        // tree.insert({
                        //     minX: p.x - 500,
                        //     minY: p.y - 500,
                        //     maxX: p.x + 500,
                        //     maxY: p.y + 500,
                        //     feature: marker,
                        // })

                    });
                });


                // // comment this in for now. Will get back later
                // utils.getMap().on('mousemove', L.Util.throttle(function (e) {
                // onmyMousemove(e);
                // }, 32));

                function onmyMousemove(e) {
                    console.log('im in')
                    mouseTarget = findFeature(e.latlng);
                }


                function findFeature(latlng) {
                    var point = {'x': latlng.lng, 'y': latlng.lat};
                    var features = tree.search({
                        minX: point.x,
                        minY: point.y,
                        maxX: point.x,
                        maxY: point.y
                    });

                    console.log('length features is ' + features.length)
                    if (features.length > 0) {
                        console.log('got it')
                    }
                    return features
                }

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
                    var targetScale = zoom <= zoomSwitch ? scaleRamp(zoom) : 1 / (2 * utils.getScale(event.zoom));
                    markerSprite.scale.set(targetScale);
                });
                masterMarkerRenderer.render(container);
                if (progress < delta) {
                    frame = requestAnimationFrame(animate);
                }
            }

            if (firstDraw || prevZoom !== zoom) {
                markerSprites.forEach(function (markerSprite) {
                    var targetScale = zoom <= zoomSwitch ? scaleRamp(zoom) : 1 / (2 * utils.getScale(event.zoom));
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

            // finally remove the preloader
            removePreloader();
            console.log('preloader removed')
        }
    }
};


// myMarkersContainer = masterMarkerContainer.getChildByName('myMarkers');
// myMarkersContainer.visible = false
// masterMarkerRenderer.render(masterMarkerContainer)


// [...new Set(all_geneData.map(d => d.Gene))].sort()