// The animation loop. Code in the animate() function is getting executed on every tick

var attributes;

function animate() {
    requestAnimationFrame(animate);
    render();
    // stats.update();
}

function setHightlightSphere(instanceId, isHighlighting) {
    if (instanceId == -1) return;
    var dummy = new THREE.Object3D();
    var loader = new THREE.TextureLoader();
    var props = {
        clearcoat: 1.0,
        clearcoatRoughness: 0,
        metalness: 0.065,
        roughness: 0.3,
        normalMap: loader.load('./src/flakes.png'),
        normalScale: new THREE.Vector2(0.1, 0.1),
        transmission: 0.0,
        transparent: true,
        side: THREE.BackSide,
    };
    var material = new THREE.MeshPhysicalMaterial(props);
    material.opacity = 1.0;
    material.normalMap.wrapS = material.normalMap.wrapT = THREE.RepeatWrapping;
    material.normalMap.repeat = new THREE.Vector2(30, 30);

    var highlighter = new THREE.InstancedMesh(
        //provide geometry
        new THREE.SphereBufferGeometry(1, 36, 18),

        //provide material
        material,

        //how many instances to allocate
        1
    );

    var coords = NON_ZERO_CELLS[instanceId].sphere_position,
        scales = NON_ZERO_CELLS[instanceId].sphere_scale,
        rot = NON_ZERO_CELLS[instanceId].sphere_rotation,
        // color = new THREE.Color("yellow");
        color = NON_ZERO_CELLS[instanceId].color;
    dummy.position.set(coords.x, coords.y, coords.z);
    dummy.scale.set(scales.x*1.02, scales.y*1.02, scales.z*1.02);
    dummy.rotation.set(rot.x, rot.y, rot.z);
    dummy.updateMatrix();
    highlighter.name = 'cell_highlight';
    highlighter.setMatrixAt(0, dummy.matrix);
    highlighter.setColorAt(0, new THREE.Color(color.r, color.g, color.b));
    highlighter.receiveShadow = false;
    highlighter.castShadow = true;

    if (isHighlighting) {
        SCENE.add(highlighter)
    }
    // else {
    //     SCENE.remove(highlighter)
    // }
    // // instancedMesh.geometry.setPositionAt(i, trsCache[i].position);
    // // instancedMesh.geometry.setScaleAt(i, uScale ? ss : trsCache[i].scale);
}

function add_highlight_sphere(instanceId) {
    if (instanceId != PREV_INSTANCE_ID) {
        remove_highlight_sphere();
        setHightlightSphere(instanceId, true);

        $('html,body').css('cursor', 'pointer');
    }
}

function remove_highlight_sphere() {
    // 1. restore previous id
    PREV_INSTANCE_ID = -1;

    // 2. remove the mesh from the scene
    SCENE.children
        .filter(d => d.name === 'cell_highlight')
        .forEach(d => SCENE.remove(d))

    // restore the mouse cursor
    $('html,body').css('cursor', 'default');
}



function render() {
    RAYCASTER.setFromCamera(MOUSE, CAMERA);

    var target = SCENE.children.filter(d => (d.type === 'Points') & (d.visible));
    const intersects = RAYCASTER.intersectObjects(target, true);

    if (intersects.length > 0) {
        var particles = SCENE.children.filter(d => (d.type === 'Points') & (d.visible) & (d.name === intersects[0].object.name))[0]
        if (particles) {
            $('html,body').css('cursor', 'pointer');
            var geometry = particles.geometry;
            attributes = geometry.attributes;


            if (INTERSECTED.index !== intersects[0].index && INTERSECTED.name !== intersects[0].object.name) {
                if (INTERSECTED.uuid) {
                    var previous_particles = SCENE.children.filter(d => (d.type === 'Points') & (d.uuid === INTERSECTED.uuid))[0];
                    previous_particles.geometry.attributes.size.array[INTERSECTED.index] = CONFIGSETTINGS.particle_size;
                    previous_particles.geometry.attributes.size.needsUpdate = true;

                    previous_particles.geometry.attributes.mycolor.array[3 * INTERSECTED.index] = INTERSECTED.rgb[0];
                    previous_particles.geometry.attributes.mycolor.array[3 * INTERSECTED.index + 1] = INTERSECTED.rgb[1];
                    previous_particles.geometry.attributes.mycolor.array[3 * INTERSECTED.index + 2] = INTERSECTED.rgb[2];
                    previous_particles.geometry.attributes.mycolor.needsUpdate = true;

                    // previous_particles.geometry.attributes.position.array[3 * INTERSECTED.index] = INTERSECTED.position[0];
                    // previous_particles.geometry.attributes.position.array[3 * INTERSECTED.index + 1] = INTERSECTED.position[1];
                    // previous_particles.geometry.attributes.position.array[3 * INTERSECTED.index + 2] = INTERSECTED.position[2];
                    // previous_particles.geometry.attributes.position.needsUpdate = true;
                }

                attributes.size.array[INTERSECTED.index] = CONFIGSETTINGS.particle_size;
                INTERSECTED.rgb = [attributes.mycolor.array[3 * intersects[0].index],
                    attributes.mycolor.array[3 * intersects[0].index + 1],
                    attributes.mycolor.array[3 * intersects[0].index + 2],
                ];
                attributes.mycolor.array[3 * intersects[0].index] = 255;
                attributes.mycolor.array[3 * intersects[0].index + 1] = 255;
                attributes.mycolor.array[3 * intersects[0].index + 2] = 0;
                attributes.mycolor.needsUpdate = true;

                INTERSECTED.position = [attributes.position.array[3 * intersects[0].index],
                    attributes.position.array[3 * intersects[0].index + 1],
                    attributes.position.array[3 * intersects[0].index + 2],
                ];

                INTERSECTED.index = intersects[0].index;
                INTERSECTED.name = intersects[0].object.name;
                INTERSECTED.uuid = intersects[0].object.uuid;

                attributes.size.array[INTERSECTED.index] = CONFIGSETTINGS.particle_size * 1.25;
                attributes.size.needsUpdate = true;

                TOOLTIP.visible = true;
                $('.pointlabel').text(INTERSECTED.name);
                fs = font_ramp(intersects[0].distance);
                m = margin(intersects[0].distance);
                $('.pointlabel').css("font-size", fs + "px");

                $('.pointlabel').css("margin-top", "0px");
                $('.pointlabel').css("margin-left", "0px");
                $('.pointlabel').css("margin-top", -m + "px");
                $('.pointlabel').css("margin-left", m + "px");
                console.log('Particle depth is: ' + intersects[0].distance);
                console.log('Font size is: ' + font_ramp(intersects[0].distance));
                TOOLTIP.position.x = attributes.position.array[3 * intersects[0].index];
                TOOLTIP.position.y = attributes.position.array[3 * intersects[0].index + 1];
                TOOLTIP.position.z = attributes.position.array[3 * intersects[0].index + 2];

            }

        }
    } else if (INTERSECTED.index !== null && attributes) {
        $('html,body').css('cursor', 'default');
        $('.pointlabel').text('')
        TOOLTIP.visible = false;

        attributes.size.array[INTERSECTED.index] = CONFIGSETTINGS.particle_size;
        attributes.size.needsUpdate = true;

        // attributes.position.array[3 * INTERSECTED.index] = INTERSECTED.position[0];
        // attributes.position.array[3 * INTERSECTED.index + 1] = INTERSECTED.position[1];
        // attributes.position.array[3 * INTERSECTED.index + 2] = INTERSECTED.position[2];
        // attributes.position.needsUpdate = true;

        attributes.mycolor.array[3 * INTERSECTED.index] = INTERSECTED.rgb[0];
        attributes.mycolor.array[3 * INTERSECTED.index + 1] = INTERSECTED.rgb[1];
        attributes.mycolor.array[3 * INTERSECTED.index + 2] = INTERSECTED.rgb[2];
        attributes.mycolor.needsUpdate = true;

        INTERSECTED.index = null;
        INTERSECTED.name = null;
        INTERSECTED.uuid = null;
        INTERSECTED.rgb = null;
    }

    const intersection = RAYCASTER.intersectObject(INSTANCEDMESH.front_face.instancedMesh);
    if (intersection.length > 0) {
        var instanceId = intersection[0].instanceId;
        var cell_label = NON_ZERO_CELLS[instanceId].Cell_Num;
        console.log('Hovering over cell: ' + cell_label)

        INSTANCEDMESH.front_face.instancedMesh.visible = false;
        SCENE.children.filter(d => (d.type === 'Points') & (d.name !== 'glyph_highlighting')).forEach(d => d.visible=false);
        var _target_spots = ALL_GENEDATA.filter(d => d.neighbour===cell_label);
        var target_spots = groupBy(_target_spots, 'Gene');
        if (CTRL_KEY_PRESSED) {
            if (cell_label !== PREVIOUS_CELL_LABEL){
                console.log(target_spots);
                var target_genes = Object.keys(target_spots);
                var temp_arr = [];
                for (var i = 0; i < target_genes.length; i++) {
                    var g = target_genes[i],
                        temp = new Float32Array(target_spots[g].map(d => centralise_coords([d.x, d.y, d.z], CONFIGSETTINGS)).flat());
                    temp_arr.push(temp)
                }
                HIGHLIGHTING_POINTS = target_genes.map((d, i) => my_particles(temp_arr[i], d));
                HIGHLIGHTING_POINTS.forEach(d => d.name = 'glyph_highlighting');
                HIGHLIGHTING_POINTS.map(d => SCENE.add(d));
                PREVIOUS_CELL_LABEL = cell_label;

                // make lines
                // SCENE.children.filter(d => d.type === "Line").forEach(el => SCENE.remove(el))
                // var centroid = [NON_ZERO_CELLS[instanceId].X, NON_ZERO_CELLS[instanceId].Y, NON_ZERO_CELLS[instanceId].Z]
                // var out = make_line(target_spots, centroid);
                // out.map(d => SCENE.add(d))
            }
            else{
                SCENE.children.filter(d => d.name === 'glyph_highlighting').forEach(d => d.visible=true);

            }

            // var points = GENEPANEL.map((d, i) => my_particles(spots[i], d));
            // points.map(d => SCENE.add(d));
        }
        else{
            INSTANCEDMESH.front_face.instancedMesh.visible = true;
            SCENE.children.filter(d => d.type === 'Points').forEach(d => d.visible=true)
            SCENE.children.filter(d => d.name === 'glyph_highlighting').forEach(d => SCENE.remove(d));
            HIGHLIGHTING_POINTS = null
            PREVIOUS_CELL_LABEL = null

            donutchart(NON_ZERO_CELLS[instanceId]);
            renderDataTable(NON_ZERO_CELLS[instanceId])

            // make lines
            SCENE.children.filter(d => d.type === "Line").forEach(el => SCENE.remove(el))
            var centroid = [NON_ZERO_CELLS[instanceId].X, NON_ZERO_CELLS[instanceId].Y, NON_ZERO_CELLS[instanceId].Z]
            var out = make_line(target_spots, centroid);
            out.map(d => SCENE.add(d))

        }

        // INSTANCEDMESH.back_face.instancedMesh.material.opacity = 1.0;
        add_highlight_sphere(instanceId);
        // add_highlight_glyphs(instanceId);
        PREV_INSTANCE_ID = instanceId;
    } else {
        remove_highlight_sphere();
        // remove_highlight_glyphs();
    }

    RENDERER.render(SCENE, CAMERA);
    LABEL_RENDERER.render(SCENE, CAMERA)
}

function make_line(obj, centroid){
    var arr = Object.entries(obj).map(d => d[1]).flat()
    var out = arr.map(d => {
        return make_line_helper(d, centroid)
    });
    return out
}

function remove_line(){
    SCENE.children.filter(d => d.type === "Line").forEach(el => SCENE.remove(el))
}

function make_line_helper(d, centroid) {
    var centre = centralise_coords([centroid[0], centroid[1], centroid[2]], CONFIGSETTINGS).flat();
    var xyz = centralise_coords([d.x, d.y, d.z], CONFIGSETTINGS).flat();
    var points = [];
    points.push(
        new THREE.Vector3(...xyz),
        new THREE.Vector3(...centre),
    )
    var geometry = new THREE.BufferGeometry().setFromPoints(points);
    // CREATE THE LINE
    var line = new THREE.Line(
        geometry,
        new THREE.LineBasicMaterial({
            // color: getColor(d.Gene)
            color: GLYPH_MAP.get(d.Gene).color
        })
    );
    return line
}

function font_ramp(z){
    return z > 1000? 10:
        z > 500? 15:
            z > 250? 20:
                z > 125? 25:
                    z > 65? 30: 35;
    // var range = (CAMERA.far - z) / (CAMERA.far - CAMERA.near);
}

function margin(depth){
    var k;
    if (depth < 280){
        k = 1.75
    }
    if (depth < 180){
        k = 1.5
    }
    else if (depth < 120)
    {
        k = 1.25
    }
    else {
        k = 2.0
    }
    var x = 8000/depth * k;
    return Math.ceil(x/2)
}