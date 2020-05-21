function drawCellSprites(resources, markers, map) {
    var textures = [resources.plane.texture, resources.circle.texture, resources.bicycle.texture];
    var focusTextures = [resources.focusPlane.texture, resources.focusCircle.texture, resources.focusBicycle.texture];

    var pixiLayer = cellOverlay(markers, textures, focusTextures);
    // pixiLayer().addTo(map);
    return pixiLayer()
}

function cellOverlay(markers, textures, focusTextures) {
    return function () {
        var pixiContainer = new PIXI.Container();
        var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;

        var _drawCallback = drawCallback(markers, textures, focusTextures);
        return L.pixiOverlay(_drawCallback, pixiContainer, {
            doubleBuffering: doubleBuffering,
            destroyInteractionManager: true
        });
    };


    function drawCallback(markers, textures, focusTextures) {
        var firstDraw = true;
        var prevZoom;
        var markerSprites = [];
        var colorScale = d3.scaleLinear()
            .domain([0, 50, 100])
            .range(["#c6233c", "#ffd300", "#008000"]);

        var frame = null;
        var focus = null;
        return function (utils, event) {
            var zoom = utils.getMap().getZoom();
            if (frame) {
                cancelAnimationFrame(frame);
                frame = null;
            }
            var container = utils.getContainer();
            var renderer = utils.getRenderer();
            var project = utils.latLngToLayerPoint;
            var scale = utils.getScale();
            var invScale = 1 / scale;
            var initialScale = invScale / 4;
            if (firstDraw) {
                prevZoom = zoom;
                markers.forEach(marker => {
                    var tempSpot = [marker.Y, marker.X],
                        lp = dapiConfig.t.transform(L.point(tempSpot)),
                        coords = project([lp.x, lp.y]);
                    markerSprite = createMarker(coords, textures, scale);
                    var tint = d3.color(colorScale(marker.avancement || Math.random() * 100)).rgb();
                    markerSprite.tint = 256 * (tint.r * 256 + tint.g) + tint.b;
                    container.addChild(markerSprite);
                    markerSprites.push(markerSprite);
                    markerSprite.legend = marker.city || marker.label;
                });

                // create zoom quad tree
                let zoomQuadTree = {};
                for (var z = map.getMinZoom(); z <= map.getMaxZoom(); z++) {
                    const initialRadius = 8 / utils.getScale(z);
                    zoomQuadTree[z] = solveCollision(markerSprites, {r0: initialRadius, zoom: z});
                }

                var self = this;
                map.on('mousemove', L.Util.throttle(function (e) {
                    onMouseMove(e);
                    // console.log('mouse move fired')
                }, 32));

                function onMouseMove(e) {
                    // console.log('clicked');
                    var redraw = false;
                    var my_target = findMarker(project(e.latlng), zoomQuadTree, zoom);
                    if (my_target) {
                        console.log('Event fired, marker found');
                        my_target.texture = textures[1]
                        redraw = true
                    }
                    if (my_target) {
                        L.DomUtil.addClass(self._container, 'leaflet-interactive');
                    } else {
                        L.DomUtil.removeClass(self._container, 'leaflet-interactive');
                    }
                    // console.log(my_target);
                    if (redraw) utils.getRenderer().render(container)
                }

            } // first draw ends here


            if (firstDraw || prevZoom !== zoom) {
                markerSprites.forEach(markerSprite => {
                    rescaleMarker(markerSprite, invScale, zoom)
                });
            }

            firstDraw = false;
            prevZoom = zoom;
            renderer.render(container);
        }
    }

    function createMarker(coords, textures, scale) {
        var index = Math.floor(Math.random() * textures.length);
        var sprite = new PIXI.Sprite(textures[index]);
        sprite.textureIndex = index;
        sprite.x0 = coords.x;
        sprite.y0 = coords.y;
        sprite.anchor.set(0.5, 0.5);
        sprite.scale.set(scale);
        sprite.currentScale = scale;
        return sprite;
    }

    function rescaleMarker(marker, scale, zoom, redraw) {
        const position = marker.cache[zoom];
        if (!redraw) { // 1st draw
            marker.x = position.x;
            marker.y = position.y;
            marker.scale.set((position.r * scale < 8) ? position.r / 8 : scale/16); // 16
        } else {
            marker.currentX = marker.x;
            marker.currentY = marker.y;
            marker.targetX = position.x;
            marker.targetY = position.y;
            marker.currentScale = marker.scale.x;
            marker.targetScale = (position.r * scale < 8) ? position.r / 16 : scale; // 16
        }
    }

    function findMarker(layerPoint, quad, zoom) {
        const quadTree = quad[zoom];
        const maxR = quadTree.rMax;
        var marker;
        var found = false;
        quadTree.visit((quad, x1, y1, x2, y2) => {
            if (!quad.length) {
                const dx = quad.data.x - layerPoint.x;
                const dy = quad.data.y - layerPoint.y;
                const r = quad.data.scale.x * 8; // 16;
                if (dx * dx + dy * dy <= r * r) {
                    marker = quad.data;
                    found = true;
                }
            }
            return found ||
                x1 > layerPoint.x + maxR ||
                x2 + maxR < layerPoint.x ||
                y1 > layerPoint.y + maxR ||
                y2 + maxR < layerPoint.y;
        });
        return marker;
    }

}
