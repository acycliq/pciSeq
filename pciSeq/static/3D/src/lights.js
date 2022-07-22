function iniLights() {
    // LIGHTS
    SCENE.add(new THREE.AmbientLight(0x666666));

    var light = new THREE.DirectionalLight(0xfcfcfc, 1.0);
    light.position.set(90, 120, 5);
    light.position.multiplyScalar(1.3);
    light.castShadow = true;
    // light.shadowCameraVisible = true;
    var d = 200;
    light.shadow.camera.left = -0.2 * d;
    light.shadow.camera.right = 0.2 * d;
    light.shadow.camera.top = d;
    light.shadow.camera.bottom = -d;
    light.shadow.camera.far = 400;
    // light.shadow.bias = -0.01; // You may need to tweak this to avoid artifacts if the mesh is receiving shadows
    // light.shadowDarkness = 0.2;
    SCENE.add(light);

    // var helper = new THREE.CameraHelper(light.shadow.camera);
    // scene.add(helper);


    var light_2;
    light_2 = new THREE.DirectionalLight(0xfcfcfc, 1.0);
    light_2.position.set(-90, -120, -5);
    light_2.position.multiplyScalar(1.3);
    light_2.castShadow = true;
    // light.shadowCameraVisible = true;
    var d = 200;
    light_2.shadow.camera.left = -0.25 * d;
    light_2.shadow.camera.right = 0.25 * d;
    light_2.shadow.camera.top = d;
    light_2.shadow.camera.bottom = -d;
    light_2.shadow.camera.far = 400;
    // light_2.shadow.bias = -0.01; // You may need to tweak this to avoid artifacts if the mesh is receiving shadows
    // light.shadowDarkness = 0.2;
    SCENE.add(light_2);

    // var helper_2 = new THREE.CameraHelper(light_2.shadow.camera);
    // scene.add(helper_2);

    add_envmap();

    // spotLight = new THREE.SpotLight(0xffffff);
    // spotLight.name = 'Spot Light';
    // spotLight.angle = Math.PI / 5;
    // spotLight.penumbra = 0.3;
    // spotLight.position.set(10, 60, -40);
    // spotLight.castShadow = true;
    // spotLight.shadow.camera.near = 80;
    // spotLight.shadow.camera.far = 150;
    // spotLight.shadow.camera.left = -d;
    // spotLight.shadow.camera.right = d;
    // spotLight.shadow.camera.top = d;
    // spotLight.shadow.camera.bottom = -d;
    // spotLight.shadow.mapSize.width = 512;
    // spotLight.shadow.mapSize.height = 512;
    // scene.add(spotLight);
    //
    // var spotLight_helper = new THREE.CameraHelper(spotLight.shadow.camera);
    // scene.add(spotLight_helper);

}

function add_envmap() {
    new THREE.RGBELoader()
        .setDataType(THREE.UnsignedByteType)
        .setPath('https://threejs.org/examples/textures/equirectangular/')
        .load('royal_esplanade_1k.hdr', function (texture) {

            var envMap = pmremGenerator.fromEquirectangular(texture).texture;

            // scene.background = envMap;
            SCENE.environment = envMap;

            texture.dispose();
            pmremGenerator.dispose();
        })
    var pmremGenerator = new THREE.PMREMGenerator(RENDERER);
    pmremGenerator.compileEquirectangularShader();

}
