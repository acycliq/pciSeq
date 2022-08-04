// Put here event listeners and other stuff to be called when everything else on the scene is ready

function postScene() {

    // adjust the scene if the browser's window is resized
    window.addEventListener('resize', onWindowResize);

    // mouse move
    window.addEventListener('mousemove', onMouseMove, false);

    document.addEventListener('keydown', (event) => {
        console.log(event); // all event related info
        console.log(event.type);
        console.log(event.key);
        console.log(event.code);
        if (event.ctrlKey){
            CTRL_KEY_PRESSED = true
        }
    });

    document.addEventListener('keyup', (event) => {
        console.log(event); // all event related info
        console.log(event.type);
        console.log(event.key);
        console.log(event.code);
        if (!event.ctrlKey){
            console.log('CTRL UP')
            CTRL_KEY_PRESSED = false
        }
    });

    // mouse wheel
    // window.addEventListener('mousewheel', MouseWheelHandler, false);

    // finally remove the preloader
    removePreloader();

    // and show the gene panel button
    legendControl();

    //
    populate_jstree()

    function onMouseMove(event) {
        // calculate mouse position in normalized device coordinates
        // (-1 to +1) for both components

        MOUSE.x = (event.clientX / window.innerWidth) * 2 - 1;
        MOUSE.y = -(event.clientY / window.innerHeight) * 2 + 1;
    }

    function onWindowResize() {
        // Update camera
        CAMERA.aspect = window.innerWidth / window.innerHeight;
        CAMERA.updateProjectionMatrix();

        // Update renderer
        RENDERER.setSize(window.innerWidth, window.innerHeight);
        RENDERER.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    }

    function MouseWheelHandler(e) {
        // this will zoom to the mouse cursor (but it breaks panning!!)
        console.log('mouse wheel event detected')
        var vector = new THREE.Vector3(MOUSE.x, MOUSE.y, 1 );
        vector.unproject(CAMERA);
        vector.sub(CAMERA.position);
        var factor = 1000.0;
        CAMERA.position.addVectors(CAMERA.position, vector.setLength(factor));
        CONTROLS.target.addVectors(CONTROLS.target, vector.setLength(factor));
    }

    function populate_jstree(){
        // var cellClasses = [...new Set(CELL_DATA.map(d => d.topClass))].sort()
        $('#jstree_id').jstree(true).settings.core.data = [tree(TOPCLASSES)]
        $('#jstree_id').jstree(true).refresh();

        $("#jstree_id").jstree("close_all");
    }


    // var paramsGUI = {
    //     numSpots: 0,
    //     numCells: 0,
    // };
    //
    // gui = new dat.GUI();
    //
    // var numSpots = [0, 100000, 1000000, 5000000, 10000000, 20000000],
    //     numCells = [0, 100, 1000, 100000, 150000, 200000];
    // gui.add(paramsGUI, 'numSpots', numSpots).name('Num simulated spots').onChange(onSelectCounts);
    // gui.add(paramsGUI, 'numCells', numCells).name('Num simulated cells').onChange(onSelectCounts);
    //
    // function onSelectCounts() {
    //     console.log('Selected: ' + paramsGUI.numSpots + ' number of spots');
    //     var spots_xyz,
    //         cells_xyz;
    //     if (+paramsGUI.numSpots || +paramsGUI.numCells) {
    //         +paramsGUI.numSpots ? spots_xyz = simulate_spots(+paramsGUI.numSpots) : spots_xyz = SPOTS_ARR;
    //         +paramsGUI.numCells ? cells_xyz = get_sim_cell_xyz(+paramsGUI.numCells) : cells_xyz = CELLS_ARR;
    //         // console.log(spots_xyz);
    //     } else {
    //         spots_xyz = SPOTS_ARR;
    //         cells_xyz = CELLS_ARR;
    //     }
    //
    //     iniScene();
    //     iniLights();
    //     iniContent(spots_xyz, cells_xyz);
    // }
    //
    // gui.open();
}