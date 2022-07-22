function iniScene() {

    // Canvas
    CANVAS = document.querySelector('canvas.webgl')

    var container = document.createElement('div');
    document.body.appendChild(container);

    var near = 20,
        far = 80000;
    CAMERA = new THREE.PerspectiveCamera(25, window.innerWidth / window.innerHeight, near, far);
    CAMERA.position.set(24000, 1665, 13500);

    SCENE = new THREE.Scene();
    // SCENE.background = new THREE.Color(0xdddddd);

    // ground
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

    // Renderer
    RENDERER = new THREE.WebGLRenderer({
        canvas: CANVAS,
        antialias: true,
    });
    RENDERER.setSize(window.innerWidth, window.innerHeight);
    RENDERER.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    RENDERER.shadowMap.enabled = true;
    RENDERER.shadowMap.type = THREE.BasicShadowMap;

    const tooltipDiv = document.createElement('div');
    tooltipDiv.className = 'pointlabel';
    tooltipDiv.textContent = 'Tooltip';
    tooltipDiv.style.marginTop = '-1em';
    TOOLTIP = new THREE.CSS2DObject(tooltipDiv);
    TOOLTIP.position.set(0, 0, 0);
    SCENE.add(TOOLTIP);

    LABEL_RENDERER = new THREE.CSS2DRenderer();
    LABEL_RENDERER.setSize(window.innerWidth, window.innerHeight);
    LABEL_RENDERER.domElement.style.position = 'absolute';
    LABEL_RENDERER.domElement.style.top = '0px';
    document.body.appendChild(LABEL_RENDERER.domElement);

    // Controls
    CONTROLS = new THREE.OrbitControls(CAMERA, LABEL_RENDERER.domElement);
    // controls.enableDamping = true;

    // stats = new Stats();
    // container.appendChild(stats.dom);

    // var axes = createAxes(1000, SCENE);
    // SCENE.add(axes)

}


