function make_cells(cellData) {
    var front_props = {
            side: THREE.FrontSide,
            opacity: 0.4,
            name: 'front_mesh'
        },
        back_props = {
            side: THREE.BackSide,
            opacity: 0.9,
            name: 'back_mesh'
        };

    cellData = cellData.filter(d => d.color.r + d.color.g + d.color.b !== 0);
    var front_face = ellipsoids(cellData, front_props),
        back_face = ellipsoids(cellData, back_props);
        cells = {};
    cells.front_face = front_face;
    cells.back_face = back_face;
    return cells
}

function ellipsoids(data, props) {
    var counts = data.length,
        loader = new THREE.TextureLoader();

    var img_width = CONFIGSETTINGS.img_width,
        img_height = CONFIGSETTINGS.img_height,
        img_depth = CONFIGSETTINGS.img_depth;

    const flakesTexture = loader.load('./src/flakes.png')
    const base_props = {
        clearcoat: 1.0,
        clearcoatRoughness: 0,
        metalness: 0.065,
        roughness: 0.3,
        normalMap: flakesTexture,
        normalScale: new THREE.Vector2(0.3, 0.3),
        transmission: 0.0,
        transparent: true,
        // envMap: true,
    };
    var material = new THREE.MeshPhysicalMaterial(base_props);

    material.side = props.side;
    material.color = props.color;
    material.opacity = props.opacity;
    material.normalMap.wrapS = material.normalMap.wrapT = THREE.RepeatWrapping;
    material.normalMap.repeat = new THREE.Vector2(30, 30)


    var uScale = 0;
    var widthSegments = 8,
        heightSegments = 0.5 * widthSegments;
    var geometry =  new THREE.SphereBufferGeometry(1, widthSegments, heightSegments);
    var _n = geometry.index.count/3;
    console.log('triangles: ' + (_n * counts).toLocaleString());
    var INSTANCEDMESH = new THREE.InstancedMesh(
        //provide geometry
        geometry,

        //provide material
        material,

        //how many instances to allocate
        counts,

        //is the scale known to be uniform, will do less shader work, improperly applying this will result in wrong shading
        !!uScale
    );

    var dummy = new THREE.Object3D();
    var bbox_items = [];
    var tree = new RBush3D.RBush3D(16)
    console.log('tic')
    var temp_obj = new THREE.Mesh(new THREE.SphereGeometry(), new THREE.MeshBasicMaterial());
    for (var i = 0; i < counts; i++) {
        var coords = data[i].sphere_position,
            scales = data[i].sphere_scale,
            rot = data[i].sphere_rotation,
            color =  data[i].color;
        dummy.position.set(coords.x, coords.y, coords.z);
        dummy.scale.set(scales.x, scales.y, scales.z);
        dummy.rotation.set(rot.x, rot.y, rot.z);
        dummy.updateMatrix();
        INSTANCEDMESH.name = props.name;
        INSTANCEDMESH.setMatrixAt(i, dummy.matrix);
        INSTANCEDMESH.setColorAt(i, new THREE.Color( color.r, color.g, color.b ));
        temp_obj.applyMatrix4(dummy.matrix)
        var bbox = new THREE.Box3().setFromObject(temp_obj);
        bbox_items.push({minX: bbox.min.x, minY: bbox.min.y, minZ: bbox.min.z, maxX: bbox.max.x, maxY: bbox.max.y, maxZ: bbox.max.z, name: 'item_' + i});
        // instancedMesh.geometry.setPositionAt(i, trsCache[i].position);
        // instancedMesh.geometry.setScaleAt(i, uScale ? ss : trsCache[i].scale);
    }
    console.log('toc')
    INSTANCEDMESH.instanceColor.needsUpdate = true;
    INSTANCEDMESH.visible = true;
    INSTANCEDMESH.castShadow = true;
    INSTANCEDMESH.receiveShadow = false; // if you turn this on, you may need to tweak the light shadow bias to avoid artifacts
    // objectWrapper.add(instancedMesh);

    tree.load(bbox_items);


    function LOD_ramp() {
        // camera.position.distanceTo(scene.position) < 300? mesh_LOD(): null
        if (camera.position.distanceTo(scene.position) < 300) {
            console.log('Less than 300')
        }

    }


    var ellipsoidData = {};
    ellipsoidData.instancedMesh = INSTANCEDMESH;
    ellipsoidData.LOD_ramp = LOD_ramp;

    return ellipsoidData
}

function count_triangles(m){
    // input m is the mesh
    var _n = m.geometry.index.count/3;
    var count = m.count;
    console.log('triangles: ' + (_n * count).toLocaleString());
}