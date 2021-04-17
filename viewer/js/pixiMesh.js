function pixiLayerMesh(topo, objName) {
		return function() {
			var firstDraw = true;
			var prevZoom;
			var pixiContainer = new PIXI.Graphics();

			var meshAlphaScale = d3.scaleLinear()
					.domain([9, 12])
					.range([0.6, 1]);
			meshAlphaScale.clamp(true);
			var doubleBuffering = /iPad|iPhone|iPod/.test(navigator.userAgent) && !window.MSStream;
			var mesh, circMesh;
			return L.pixiOverlay(function (utils) {
				var zoom = utils.getMap().getZoom();
				var container = utils.getContainer();
				var renderer = utils.getRenderer();
				var gl = renderer.gl;
				var project = utils.latLngToLayerPoint;
				var scale = utils.getScale();
				var invScale = 1 / scale;
				var self = this;
				if (firstDraw) {
					(function () {
						if (renderer.type === PIXI.RENDERER_TYPE.WEBGL) {
							gl.blendFunc(gl.ONE, gl.ZERO);
						}
						// topo.arcs.forEach(function (arc) {
						// 	arc.forEach(function (position) {
						// 		var proj = project([position[1], position[0]]);
						// 		position[0] = proj.x;
						// 		position[1] = proj.y;
						// 	});
						// });
						var interiors = topojson.mesh(topo, topo.objects[objName], function (a, b) {
							return a !== b && a.properties.ref === b.properties.ref;
						});
						var circs = topojson.mesh(topo, topo.objects[objName], function (a, b) {
							return a !== b && a.properties.cell_id !== b.properties.cell_id;
						});
						topo = null;

						prevZoom = zoom;


						if (renderer.type === PIXI.RENDERER_TYPE.WEBGL) {
							(function () {
								mesh = new PIXI.Container();
								circMesh = new PIXI.Container();

								var memo = Object.create(null);
								var newIndex = 0;
								var meshVertices = [];
								var meshIndices = [];
								var iMax, iMin;

								function meshCreate(meshVertices, meshIndices, target, color) {
									color = 0x16A085
									var partialMesh = new PIXI.mesh.Mesh(null, new Float32Array(meshVertices), null, new Uint16Array(meshIndices));
									partialMesh.tint = color;
									target.addChild(partialMesh);
								}

								function meshCb(polygon) {
									if (newIndex > 60000) {
										memo = Object.create(null);
										meshCreate(meshVertices, meshIndices, mesh, 0x333333);
										newIndex = 0;
										meshVertices = [];
										meshIndices = [];
									}
									var indices = polygon.map(function (point) {
										var key = point[0] + '#' + point[1];
										var index = memo[key];
										if (index !== undefined) return index;
										else {
											var index = memo[key] = newIndex++;
											meshVertices.push(point[0], point[1]);
											return index;
										}
									});
									iMax = polygon.length - 1;
									iMin = 0;
									meshIndices.push(indices[iMax]);
									while (iMax - iMin >= 2) {
										meshIndices.push(indices[iMax--], indices[iMin++]);
									}
									if (iMax === iMin) {
										meshIndices.push(indices[iMax], indices[iMax]);
									} else meshIndices.push(indices[iMax], indices[iMin], indices[iMin]);
								}

								function circMeshCb(triangle) {
									if (newIndex > 60000) {
										memo = Object.create(null);
										meshCreate(meshVertices, meshIndices, circMesh, 0);
										newIndex = 0;
										meshVertices = [];
										meshIndices = [];
									}
									var indices = triangle.map(function (point) {
										var key = point[0] + '#' + point[1];
										var index = memo[key];
										if (index !== undefined) return index;
										else {
											var index = memo[key] = newIndex++;
											meshVertices.push(point[0], point[1]);
											return index;
										}
									});
									iMax = triangle.length - 1;
									iMin = 0;
									meshIndices.push(indices[iMax]);
									while (iMax - iMin >= 2) {
										meshIndices.push(indices[iMax--], indices[iMin++]);
									}
									if (iMax === iMin) {
										meshIndices.push(indices[iMax], indices[iMax]);
									} else meshIndices.push(indices[iMax], indices[iMin], indices[iMin]);
								}

								var point2index = {};
								var vertices = [];
								var edges = [];
								interiors.coordinates.forEach(function (arc) {
									arc.forEach(function (point, index) {
										var key = point[0] + '#' + point[1];
										var indexTo;
										if (!(key in point2index)) {
											indexTo = point2index[key] = vertices.length;
											vertices.push(point);
										} else {
											indexTo = point2index[key];
										}
										if (index > 0) {
											var prevPoint = arc[index - 1];
											var indexFrom = point2index[prevPoint[0] + '#' + prevPoint[1]];
											if (indexFrom !== indexTo) edges.push([indexTo, indexFrom]);
										}
									})
								});
								graphDraw({vertices: vertices, edges: edges}, 2 / utils.getScale(12), meshCb, Math.PI);
								meshCreate(meshVertices, meshIndices, mesh, 0xff0000);
								memo = Object.create(null);
								newIndex = 0;
								meshVertices = [];
								meshIndices = [];

								var point2index2 = {};
								var vertices2 = [];
								var edges2 = [];
								circs.coordinates.forEach(function (arc) {
									arc.forEach(function (point, index) {
										var key = point[0] + '#' + point[1];
										var indexTo;
										if (!(key in point2index2)) {
											indexTo = point2index2[key] = vertices2.length;
											vertices2.push(point);
										} else {
											indexTo = point2index2[key];
										}
										if (index > 0) {
											var prevPoint = arc[index - 1];
											var indexFrom = point2index2[prevPoint[0] + '#' + prevPoint[1]];
											edges2.push([indexTo, indexFrom]);
										}
									})
								});
								graphDraw({
									vertices: vertices2,
									edges: edges2
								}, 6 / utils.getScale(12), circMeshCb, Math.PI);
								meshCreate(meshVertices, meshIndices, circMesh, 0);
							})();
						} else {
							mesh = new PIXI.Graphics();
							mesh.lineStyle(2 / utils.getScale(12), 0x333333, 1);
							interiors.coordinates.forEach(function (path) {
								path.forEach(function (point, index) {
									if (index === 0) mesh.moveTo(point[0], point[1]);
									else mesh.lineTo(point[0], point[1]);
								});
							});
							circMesh = new PIXI.Graphics();
							circMesh.lineStyle(6 / utils.getScale(12), 0x000000, 1);
							circs.coordinates.forEach(function (path) {
								path.forEach(function (point, index) {
									if (index === 0) circMesh.moveTo(point[0], point[1]);
									else circMesh.lineTo(point[0], point[1]);
								});
							});
						}
						interiors = null;
						circs = null;

						container.addChild(mesh);
						container.addChild(circMesh);

					})();
				}
				firstDraw = false;
				mesh.visible = (zoom >= 9);
				mesh.alpha = meshAlphaScale(zoom);
				circMesh.alpha = meshAlphaScale(zoom);
				prevZoom = zoom;
				renderer.render(container);
			}, pixiContainer, {
				doubleBuffering: doubleBuffering,
				destroyInteractionManager: true
			});
		}
	};