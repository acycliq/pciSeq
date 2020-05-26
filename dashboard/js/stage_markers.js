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
            c.tint = markerTint(d)
            geneContainer_array.push(c)
        });
    }

    function markerTint(geneName) {
        var hex = glyphColor(glyphSettings().filter(d => d.gene === geneName)[0].taxonomy)
        var out = myUtils().string2hex(hex);
        return out
    }


    function myCallback(markers, textures, focusTextures, all_geneData) {
        var firstDraw = true;
        var prevZoom;
        var markerSprites = [];

        var frame = null;
        var focus = null;
        return function (utils) {
            var zoom = utils.getMap().getZoom();
            var container = utils.getContainer(); // <----- I think this is referencing the same obj as '''masterMarkerContainer'''
            masterMarkerRenderer = utils.getRenderer();
            var project = utils.latLngToLayerPoint;
            var scale = utils.getScale();
            var invScale = 1 / scale;

            var my_gene = 'Npy'
            // var pixiParticleContainer = new PIXI.ParticleContainer(all_geneData.length, {vertices: true});
            // var pixiParticleContainer = geneContainer_array.filter(d => d.name === my_gene)[0];
            // var colorScale = d3.scaleLinear()
            //     .domain([0, 50, 100])
            //     .range(["#c6233c", "#ffd300", "#008000"]);
            // var tint = d3.color(colorScale(Math.random() * 100)).rgb();
            // console.log('tint is ' + tint)
            // pixiParticleContainer.tint = 256 * (tint.r * 256 + tint.g) + tint.b;
            // pixiParticleContainer.name = 'myMarkers'
            // container.addChild(pixiParticleContainer)

            geneContainer_array.forEach(d => container.addChild(d));
            var all_genes = geneContainer_array.map(d => d.name).sort()
            if (firstDraw) {
                prevZoom = zoom;
                all_genes.forEach(gene => {
                    var pixiParticleContainer = geneContainer_array.filter(d => d.name === gene)[0];
                    var tempData = all_geneData.filter(d => d.Gene === gene)
                    tempData.forEach(function (marker) {
                        var p = dapiConfig.t.transform(L.point([marker.x, marker.y]));
                        var coords = project([p.y, p.x]);
                        var markerSprite = new PIXI.Sprite(textures[0]);
                        markerSprite.x = coords.x;
                        markerSprite.y = coords.y;
                        markerSprite.anchor.set(0.5, 0.5);
                        pixiParticleContainer.addChild(markerSprite); //<===== HERE THE MARKER IS ADDED TO THE PARTICLECONTAINER
                        markerSprites.push(markerSprite);
                        markerSprite.legend = marker.city || marker.label;
                    });

                })

                // all_geneData.forEach(function (marker) {
                //     if (marker.Gene === my_gene) {
                //         var p = dapiConfig.t.transform(L.point([marker.x, marker.y]));
                //         var coords = project([p.y, p.x]);
                //         var markerSprite = new PIXI.Sprite(textures[0]);
                //         markerSprite.x = coords.x;
                //         markerSprite.y = coords.y;
                //         markerSprite.anchor.set(0.5, 0.5);
                //         pixiParticleContainer.addChild(markerSprite); //<===== HERE THE MARKER IS ADDED TO THE PARTICLECONTAINER
                //         markerSprites.push(markerSprite);
                //         markerSprite.legend = marker.city || marker.label;
                //     }
                // });
            }
            if (firstDraw || prevZoom !== zoom) {
                markerSprites.forEach(function (markerSprite) {
                    markerSprite.scale.set(invScale / 2);
                });
            }
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