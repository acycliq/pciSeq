function render_scene(points, cells, sphere) {
    // Canvas
    const canvas = document.querySelector('canvas.webgl')

    // Scene
    scene = new THREE.Scene();
    // scene.background = new THREE.Color('white');

    // add the points to the scene
    points.map(d => scene.add(d));

    // add the cells too
    // scene.add(cells);

    scene.add(sphere)

    // const cube = new THREE.Mesh(
    //     new THREE.BoxGeometry(),
    //     new THREE.MeshBasicMaterial()
    // )
    // scene.add(cube)

    function addLight(x, y, z) {
        const color = 0xFFFFFF;
        const intensity = 1;
        const light = new THREE.DirectionalLight(color, intensity);
        light.position.set(x, y, z);
        scene.add(light);
    }

    addLight(-1, 2, 4);
    addLight(1, -1, -2);


    // Base camera
    camera = new THREE.PerspectiveCamera(125,
        window.innerWidth / window.innerHeight,
        0.0000000001, 10000);
    // camera.position.z = 3000; // vizgen
    camera.position.z = 300;  // ucl
    scene.add(camera);


    // Controls
    controls = new THREE.OrbitControls(camera, canvas);
    controls.enableDamping = true;


    // Renderer
    renderer = new THREE.WebGLRenderer({
        canvas: canvas,
    });
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));


    // Animate
    const clock = new THREE.Clock();

    const tick = () => {
        const elapsedTime = clock.getElapsedTime();

        // Update controls
        controls.update();

        // Render
        renderer.render(scene, camera);

        // Call tick again on the next frame
        window.requestAnimationFrame(tick)
    };

    tick();

    // adjust the scene if the browser's window is resized
    window.addEventListener( 'resize', onWindowResize );

    // finally remove the preloader
    removePreloader()

    // and show the gene panel button
    legendControl();
}


function onWindowResize() {
    // Update camera
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();

    // Update renderer
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
}

function hidePoints(gene){
    var points = scene.children;
    for (var i=0; i<points.length; i++){
        if (points[i].name === gene){
            points[i].visible = false;
            console.log('Gene: ' + gene + ' is switched off.')
        }
    }
//     // points.filter(d => d.name === gene)[0].visible = false;
}

function showPoints(gene){
    if (GENEPANEL.includes(gene)){
        var points = scene.children;
        points.filter(d => d.name === gene)[0].visible = true;
    }
    else{
        console.log('Gene: ' + gene + ' not in the gene panel')
    }

}

function render_scene_2() {
    function initScene() {

        camera = new THREE.PerspectiveCamera(65, window.innerWidth / window.innerHeight, 0.00001, 10000);
        camera.position.set(-915, 980, 1550);

        scene = new THREE.Scene();

        // Lights

        scene.add(new THREE.AmbientLight(0x404040));

        spotLight = new THREE.SpotLight(0xffffff);
        spotLight.name = 'Spot Light';
        spotLight.angle = Math.PI / 5;
        spotLight.penumbra = 0.3;
        spotLight.position.set(20, 20, 15);
        spotLight.castShadow = true;
        spotLight.shadow.camera.near = 8;
        spotLight.shadow.camera.far = 30;
        spotLight.shadow.mapSize.width = 1024;
        spotLight.shadow.mapSize.height = 1024;
        scene.add(spotLight);

        dirLight = new THREE.DirectionalLight(0xffffff, 1);
        dirLight.name = 'Dir. Light';
        dirLight.position.set(0, 10, 0);
        dirLight.castShadow = true;
        dirLight.shadow.camera.near = 1;
        dirLight.shadow.camera.far = 10;
        dirLight.shadow.camera.right = 15;
        dirLight.shadow.camera.left = -15;
        dirLight.shadow.camera.top = 15;
        dirLight.shadow.camera.bottom = -15;
        dirLight.shadow.mapSize.width = 1024;
        dirLight.shadow.mapSize.height = 1024;
        scene.add(dirLight);


        var geometry = new THREE.BoxGeometry(10, 0.15, 10);
        var material = new THREE.MeshPhongMaterial({
            color: 0xa0adaf,
            shininess: 150,
            specular: 0xffffff,
            shading: THREE.SmoothShading
        });

        var ground = new THREE.Mesh(geometry, material);
        ground.scale.multiplyScalar(3);
        ground.castShadow = false;
        ground.receiveShadow = true;
        // scene.add(ground);

    }

    function initMisc() {

        renderer = new THREE.WebGLRenderer();
        renderer.setClearColor(0x000000);
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.shadowMap.enabled = true;
        renderer.shadowMap.type = THREE.BasicShadowMap;

        // Mouse control
        controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.target.set(0, 2, 0);
        controls.update();

        clock = new THREE.Clock();


    }

    function initExample() {

        //setup
        options = {
            dynamic: false,
            uniformScale: false,
            material: 'phong',
            instanceNumber: 2000,
            geometry: 'sphere',
        };

        var geometries = {
            sphere: new THREE.SphereBufferGeometry(1, 12, 8),
        };

        var instanceNumbers = 100;
        var phong_material = new THREE.MeshPhongMaterial({
            color: 0xffcccc,
            shininess: 150,
            specular: 0x222222,
            shading: THREE.SmoothShading,
            // wireframe: true,
            transparent: true,
            opacity: 0.3,
        });


        //mesh wrapper
        objectWrapper = new THREE.Object3D();
        objectWrapper.position.y = 16;
        scene.add(objectWrapper);

        trsCache = [];
        console.log('Initializing object cache for 100k objects...');
        console.time('Object cache initialized.');

        for (var i = 0; i < instanceNumbers; i++) {
            trsCache.push({
                position: new THREE.Vector3(Math.random() * 90 - 90/2, -1 * Math.random()*3, Math.random() * 20 - 10).multiplyScalar(14),
                scale: new THREE.Vector3(10*(Math.random() + .5), 10*(Math.random() + .5), 10*(Math.random() + .5))
            });

        }
        console.timeEnd('Object cache initialized.');
        console.log('Initializing instanced mesh permutations...');
        console.time('Instanced mesh permutations initialized.');


        var uScale = 0;
        var instancedMesh = new THREE.InstancedMesh(
            //provide geometry
            geometries['sphere'],

            //provide material
            phong_material,

            //how many instances to allocate
            instanceNumbers,

            //is the scale known to be uniform, will do less shader work, improperly applying this will result in wrong shading
            !!uScale
        );

        var ss = new THREE.Vector3(1, 1, 1);
        var dummy = new THREE.Object3D();
        for (var i = 0; i < instanceNumbers; i++) {
            var coords = trsCache[i].position;
            var scales = trsCache[i].scale;
            dummy.position.set(coords.x, coords.y, coords.z);
            dummy.scale.set(scales.x, scales.y, scales.z);
            dummy.updateMatrix();
            instancedMesh.setMatrixAt(i, dummy.matrix);

            // instancedMesh.geometry.setPositionAt(i, trsCache[i].position);
            // instancedMesh.geometry.setScaleAt(i, uScale ? ss : trsCache[i].scale);
        }
        instancedMesh.visible = true;
        instancedMesh.castShadow = true;
        instancedMesh.receiveShadow = true;
        objectWrapper.add(instancedMesh);

        var points = GENEPANEL.map((d, i) => my_particles(coords_arr[i], d));
        // add the points to the scene
        points.map(d => scene.add(d));

        // finally remove the preloader
        removePreloader()

        // and show the gene panel button
        legendControl();
        return objectWrapper

    }


    function animate() {
        requestAnimationFrame(animate);
        renderer.render(scene, camera);

    }

    return {
        initScene: initScene,
        initMisc: initMisc,
        initExample: initExample,
        animate: animate
    }
}

