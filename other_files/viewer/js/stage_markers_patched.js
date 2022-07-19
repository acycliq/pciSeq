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
        var scale = 1 / 64;
        return z === 0 ? scaleRampHelper(z, scale) :
                z === 1 ? scaleRampHelper(z,  2 * scale) :
                    z === 2 ? scaleRampHelper(z,  2 * scale) : scaleRampHelper(z, 4 * scale)

        // return z === 0 ? 0.03 * 2**3 :
        //     z === 1 ? 0.03 * 2**3 :
        //         z === 2 ? 0.03 * 2**3 :
        //             z === 3 ? 0.03 * 2**3 :
        //                 z === 4 ? 0.03 * 2**3 :
        //                     z === 5 ?  (2**z)/7602 :
        //                         z === 6 ? 0.03 * 2**1 : // every time you zoom in, leaflet scales up by 2. Divide here by 2 to keep the marker the same as in zoom level 5
        //                             z === 7 ? 0.03 * 2**0 :
        //                                 z === 8 ? 0.03 : (2**z)/7602
    }

    function scaleRampHelper(z, scale){
        // makes a tiny dot and the its scales it up based on the map and the dapi dimensions
        // As a general remark also, keep in mind that every time you zoom in, leaflet (I think) scales up by 2.
        // Divide by 2 to keep the marker the same as size. Hence if for zoom level = 3 the  return value from
        // this function is lets say zo 10, then when to keep the same size on the screen for the dot, at zoom = 4
        // the return value should be 5
        var map_size = Math.max(...configSettings.imageSize),
            dapi_size = [configSettings.roi.x1 - configSettings.roi.x0, configSettings.roi.y1 - configSettings.roi.y0],
            max_dapi = Math.max(...dapi_size),
            c = map_size / max_dapi,
            tiny_dot = 1 / (2**z),
            dot = c * tiny_dot;
        return dot * scale
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
                var radius = 16;
                var texture = generateCircleTexture(my_color, radius, masterMarkerRenderer);
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
