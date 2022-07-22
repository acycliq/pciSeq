function load_dapi() {
    var LoadersVolume = AMI.default.Loaders.Volume;
    var HelpersStack = AMI.default.Helpers.Stack;

    // Setup loader
    var loader = new LoadersVolume(CANVAS);
    var files = ['data/dapi/dapi_image-0068.nii.gz'];

    loader.load(files)
        .then(function () {
            // merge files into clean series/stack/frame structure
            var series = loader.data[0].mergeSeries(loader.data);
            loader.free();
            loader = null;
            console.log(series[0].stack[0]);
            // be carefull that series and target stack exist!
            var stackHelper = new HelpersStack(series[0].stack[0]);
            stackHelper.border.color = 0xFFEB3B;
            gui(stackHelper);
            SCENE.add(stackHelper);
            // center camera and interactor to center of bouding box
            // for nicer experience
            var centerLPS = stackHelper.stack.worldCenter();
            CAMERA.lookAt(centerLPS.x, centerLPS.y, centerLPS.z);
            CAMERA.updateProjectionMatrix();
            CONTROLS.target.set(centerLPS.x, centerLPS.y, centerLPS.z);
        })
        .catch(function (error) {
            window.console.log('oops... something went wrong...');
            window.console.log(error);
        });

    function gui(stackHelper) {
        var stack = stackHelper.stack;
        var gui = new dat.GUI({
            autoPlace: false,
        });
        var customContainer = document.getElementById('my-gui-container');
        customContainer.appendChild(gui.domElement);
        // stack
        var stackFolder = gui.addFolder('Stack');
        // index range depends on stackHelper orientation.
        var index = stackFolder.add(
            stackHelper, 'index', 0, stack.dimensionsIJK.z - 1).step(1).listen();
    }


}