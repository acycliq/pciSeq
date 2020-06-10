function add_spots_patched(all_geneData, map) {

    function markerColor(geneName) {
        var colorCode = glyphColor(glyphSettings().filter(d => d.gene === geneName)[0].taxonomy);
        var out = myUtils().string2hex(colorCode);
        return out
    }

    function generateCircleTexture(color, radius, renderer) {
        const gfx = new PIXI.Graphics();
        const tileSize = radius * 2;
        const texture = PIXI.RenderTexture.create(tileSize, tileSize);
        gfx.beginFill(color); // color base
        gfx.alpha = 0.8;
        gfx.drawCircle(tileSize / 2, tileSize / 2, radius);
        gfx.endFill();
        renderer.render(gfx, texture);
        return texture;
    }

    function markerColor(geneName) {
        var colorCode = glyphColor(glyphSettings().filter(d => d.gene === geneName)[0].taxonomy);
        var out = myUtils().string2hex(colorCode);
        return out
    }

    const groupBy = (array, key) => {
        // from https://learnwithparam.com/blog/how-to-group-by-array-of-objects-using-a-key/
        // Return the end result
        return array.reduce((result, currentValue) => {
            // If an array already present for key, push it to the array. Else create an array and push the object
            (result[currentValue[key]] = result[currentValue[key]] || []).push(
                currentValue
            );
            // Return the current iteration `result` value, this will be taken as next iteration `result` value and accumulate
            return result;
        }, {}); // empty object is the initial value for result object
    };

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


    var markersLength = all_geneData.length;

    loader.load(function (loader, resources) {
        var textures = [resources.plane.texture, resources.circle.texture, resources.bicycle.texture];
        // var texture = textures[1];
        // var texture = generateCircleTexture(5)

        var pixiLayer = (function () {
            var zoomChangeTs = null;
            masterMarkerContainer = new PIXI.Graphics();

            // group by gene name
            var data = groupBy(all_geneData, 'Gene');

            // get all the gene names
            var geneNames = Object.keys(data).sort();

            // populate an array with empty particle containers. One for each gene
            geneNames.forEach(gene => {
                var n = data[gene].length;
                var pc = new PIXI.particles.ParticleContainer(n, {vertices: true});
                pc.anchor = {x: 0.5, y: 0.5};
                pc.x = 0;
                pc.y = 0;
                pc.name = gene;
                masterMarkerContainer.addChild(pc);
                geneContainer_array.push(pc)
            })

            // var dummy_gene = 'Sst'
            // var innerContainer = geneContainer_array.filter(d => d.name === dummy_gene)[0];

            // pixiContainer.addChild(innerContainer);
            var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;
            var initialScale;
            return L.pixiOverlay(function (utils, event) {
                var zoom = utils.getMap().getZoom();
                var container = utils.getContainer();
               masterMarkerRenderer = utils.getRenderer();
                var project = utils.latLngToLayerPoint;
                var getScale = utils.getScale;
                var invScale = 1 / getScale();


                geneNames.forEach(gene => {
                    var my_color = markerColor(gene);
                    var texture = generateCircleTexture(my_color, 16, masterMarkerRenderer)
                    var pc = geneContainer_array.filter(d => d.name === gene)[0];
                    pc.texture = texture;
                    pc.baseTexture = texture.baseTexture;
                })

                // var innerContainer = geneContainer_array.filter(d => d.name === dummy_gene)[0];
                if (event.type === 'add') {

                    initialScale = invScale / 8;
                    initialScale = 0.125;
                    var targetScale = zoom <= 7 ? scaleRamp(zoom) : 1 / (2 * utils.getScale(event.zoom));
                    for (var i = 0; i < geneNames.length; i++) {
                        var gene = geneNames[i];
                        var innerContainer = geneContainer_array.filter(d => d.name === gene)[0];
                        innerContainer.localScale = targetScale;
                        var _data = data[gene];
                        for (var j = 0; j < _data.length; j++) {
                            // our patched particleContainer accepts simple {x: ..., y: ...} objects as children:
                            var x = _data[j].x;
                            var y = _data[j].y;
                            var point = dapiConfig.t.transform(L.point([x, y]));
                            var coords = project([point.y, point.x]);
                            innerContainer.addChild({
                                x: coords.x,
                                y: coords.y,
                            });
                        }
                    }
                }

                masterMarkerRenderer.render(masterMarkerContainer);
            }, masterMarkerContainer, {
                doubleBuffering: true,
                destroyInteractionManager: true
            }); // L.pixiOverlay closes
        })();

        pixiLayer.addTo(map);

        // var ticker = new PIXI.ticker.Ticker();
        // ticker.add(function (delta) {
        //     pixiLayer.redraw({type: 'redraw', delta: delta});
        // });
        // map.on('zoomstart', function () {
        //     ticker.start();
        // });
        // map.on('zoomend', function () {
        //     ticker.stop();
        // });
        // map.on('zoomanim', pixiLayer.redraw, pixiLayer);

        removePreloader();
    });
};


// myMarkersContainer = masterMarkerContainer.getChildByName('myMarkers');
// myMarkersContainer.visible = false
// masterMarkerRenderer.render(masterMarkerContainer)


// [...new Set(all_geneData.map(d => d.Gene))].sort()