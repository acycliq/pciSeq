function add_spots_patched(all_geneData, map) {

    function markerColor(geneName) {
        var colorCode = glyphColor(glyphSettings().filter(d => d.gene === geneName)[0].taxonomy);
        var out = myUtils().string2hex(colorCode);
        return out
    }

    function generateCircleTexture(renderer, radius) {
        const gfx = new PIXI.Graphics();
        const tileSize = radius * 2;
        const texture = PIXI.RenderTexture.create(tileSize, tileSize);
        gfx.beginFill(0xFFC700); // color base
        gfx.alpha = 0.8;
        gfx.drawCircle(tileSize / 2, tileSize / 2, radius);
        gfx.endFill();
        renderer.render(gfx, texture);
        return texture;
    }

    var markersLength = all_geneData.length;

    loader.load(function (loader, resources) {
        var textures = [resources.plane.texture, resources.circle.texture, resources.bicycle.texture];
        // var texture = textures[1];
        // var texture = generateCircleTexture(5)

        var pixiLayer = (function () {
            var zoomChangeTs = null;
            var pixiContainer = new PIXI.Graphics();
            var innerContainer = new PIXI.particles.ParticleContainer(markersLength, {vertices: true, tint: true});
            // add properties for our patched particleRenderer:
            // innerContainer.texture = texture;
            // innerContainer.baseTexture = texture.baseTexture;
            innerContainer.anchor = {x: 0.5, y: 0.5};

            pixiContainer.addChild(innerContainer);
            var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;
            var initialScale;
            return L.pixiOverlay(function (utils, event) {
                var zoom = utils.getMap().getZoom();
                var container = utils.getContainer();
                var renderer = utils.getRenderer();
                var project = utils.latLngToLayerPoint;
                var getScale = utils.getScale;
                var invScale = 1 / getScale();

                var texture = generateCircleTexture(renderer, 16)
                innerContainer.texture = texture;
                innerContainer.baseTexture = texture.baseTexture;

                if (event.type === 'add') {
                    innerContainer.x = 0;
                    innerContainer.y = 0;
                    // innerContainer.tint = markerColor('Pvalb');
                    initialScale = invScale / 8;
                    innerContainer.localScale = initialScale;
                    for (var i = 0; i < markersLength; i++) {
                        // our patched particleContainer accepts simple {x: ..., y: ...} objects as children:
                        var x = all_geneData[i].x;
                        var y = all_geneData[i].y;
                        var point = dapiConfig.t.transform(L.point([x, y]));
                        var coords = project([point.y, point.x]);
                        innerContainer.addChild({
                            x: coords.x,
                            y: coords.y,
                        });
                    }
                }

                renderer.render(pixiContainer);
            }, pixiContainer, {
                doubleBuffering: true,
                destroyInteractionManager: true
            }); // L.pixiOverlay closes
        })();

        pixiLayer.addTo(map);
        removePreloader();
    });
};


// myMarkersContainer = masterMarkerContainer.getChildByName('myMarkers');
// myMarkersContainer.visible = false
// masterMarkerRenderer.render(masterMarkerContainer)


// [...new Set(all_geneData.map(d => d.Gene))].sort()