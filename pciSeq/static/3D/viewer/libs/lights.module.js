import * as THREE from "./three.js/build/three.module.js";
import {RGBELoader} from "./RGBELoader.js"

function iniLights() {
    // LIGHTS
    viewer.scene.scene.add(new THREE.AmbientLight(0x666666));

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
    viewer.scene.scene.add(light);


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
    viewer.scene.scene.add(light_2);

    // add_envmap();

}

function add_envmap() {
    new RGBELoader()
        .setDataType(THREE.UnsignedByteType)
        .setPath('https://threejs.org/examples/textures/equirectangular/')
        .load('royal_esplanade_1k.hdr', function (texture) {

            var envMap = pmremGenerator.fromEquirectangular(texture).texture;

            // scene.background = envMap;
            viewer.scene.scene.environment = envMap;

            texture.dispose();
            pmremGenerator.dispose();
        })
    // var pmremGenerator = new THREE.PMREMGenerator(RENDERER);
    // pmremGenerator.compileEquirectangularShader();

}

export default iniLights