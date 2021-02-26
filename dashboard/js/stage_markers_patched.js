function add_spots_patched(all_geneData, map) {
    // this is used to simulate leaflet zoom animation timing:
	var easing = BezierEasing(0, 0, 0.25, 1);

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
            z === 1 ? 0.150 :
                z === 2 ? 0.150 :
                    z === 3 ? 0.150 :
                        z === 4 ? 0.125 :
                            z === 5 ? 0.125/2 :
                                z === 6 ? 0.0625/2 : // every time you zoom in, leaflet scales up by 2. Divide here by 2 to keep the marker the same as in zoom level 5
                                    z === 7 ? 0.03125 : 1
    }


    var pixiLayer = (function () {
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
            geneContainer_array.push(pc) // I think I can get rid of this. I can access the pc via masterMarkerContainer.getChildByName
        });

        // pixiContainer.addChild(innerContainer);
        var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;
        var initialScale;
        var firstDraw = true;
        var prevZoom;
        return L.pixiOverlay(function (utils, event) {
            var zoom = utils.getMap().getZoom();
            var container = utils.getContainer();
            masterMarkerRenderer = utils.getRenderer();
            var project = utils.latLngToLayerPoint;
            var getScale = utils.getScale;
            var invScale = 1 / getScale();


            geneNames.forEach(gene => {
                var my_color = markerColor(gene);
                var texture = generateCircleTexture(my_color, 16, masterMarkerRenderer);
                var pc = geneContainer_array.filter(d => d.name === gene)[0];
                pc.texture = texture;
                pc.baseTexture = texture.baseTexture;
            })

            // var innerContainer = geneContainer_array.filter(d => d.name === dummy_gene)[0];
            if (firstDraw) {
                if (event.type === 'add') {

                    initialScale = invScale / 8;
                    initialScale = 0.125;
                    var targetScale = zoom <= zoomSwitch ? scaleRamp(zoom) : 1;
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
            }

            if (prevZoom !== zoom) {
                // console.log('zoom: From ' + prevZoom + ' to ' + zoom);
                geneNames.forEach(gene => {
                    var innerContainer = masterMarkerContainer.getChildByName(gene);
                    innerContainer.localScale = scaleRamp(zoom)
                })
            }
            firstDraw = false;
            prevZoom = zoom;


            masterMarkerRenderer.render(masterMarkerContainer);
        }, masterMarkerContainer, {
            doubleBuffering: true,
            destroyInteractionManager: true
        }); // L.pixiOverlay closes
    })();

    pixiLayer.addTo(map);

    // All done, hide the preloader now.
    removePreloader();

};
