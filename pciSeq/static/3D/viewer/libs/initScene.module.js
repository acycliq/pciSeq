import * as THREE from "./three.js/build/three.module.js";
import make_cells_2 from "./stage_cells.module.js";
import {tree, myjsTree} from "./stage_cells.module.js";
import iniLights from "./lights.module.js";

var last_visited = 0
function initScene(cellData){
    cells = make_cells_2(cellData)
    viewer.scene.scene.add(cells.front_face.instancedMesh);

    iniLights()

    // viewer.renderer.domElement.addEventListener('mousemove', onMouseMove, false)
    // add a tiny throttle. The callback onMouseMove wont be called until after 5milliseconds
    // have passed.
    viewer.renderer.domElement.addEventListener('mousemove', throttle(onMouseMove, 5));


    attach_tree_control(cellData)

    // all done, remove the preloader
    removePreloader()

    clearScreen()

}

function attach_tree_control(cellData){
    var classNames = cellData.map(d => d.top_class).filter(d=>d >= 0)
    classNames = [...new Set(classNames)];
    classNames = classNames.sort()
    var treeData = tree(classNames)
    myjsTree(treeData)
}

function throttle(callback, interval) {
    // From https://programmingwithmosh.com/javascript/javascript-throttle-and-debounce-patterns/
    // similar to: https://stackoverflow.com/questions/23181243/throttling-a-mousemove-event-to-fire-no-more-than-5-times-a-second
  let enableCall = true;

  return function(...args) {
    if (!enableCall) return;

    enableCall = false;
    callback.apply(this, args);
    setTimeout(() => enableCall = true, interval);
  }
}

export function groupBy(array, key){
    // from https://learnwithparam.com/blog/how-to-group-by-array-of-objects-using-a-key/
    // Return the end result
    return array.reduce((result, currentValue) => {
        // If an array already present for key, push it to the array. Else create an array and push the object
        (result[currentValue[key]] = result[currentValue[key]] || []).push(
            currentValue
        );
        // Return the current iteration `result` value, this will be taken as next iteration `result` value and accumulate
        return result;
    }, {}); // empty object is the initial value for result object
};

function onMouseMove(event) {
    const mouse = {
        x: (event.clientX / viewer.renderer.domElement.clientWidth) * 2 - 1,
        y: -(event.clientY / viewer.renderer.domElement.clientHeight) * 2 + 1,
    }

    // console.log(mouse)
    const raycaster = new THREE.Raycaster()

    raycaster.setFromCamera(mouse, scene.getActiveCamera())

    const intersects = raycaster.intersectObjects([cells.front_face.instancedMesh])

    if (intersects.length > 0) {
        if (intersects[0].distance < 2000){
            var instanceId = cellData[intersects[0].instanceId];
            if (last_visited !== instanceId.label){
                // remove the lines from the last visited cell and draw the ones over the new cell
                // I am removing the lines twice. First all the lines (no matter the cell) and then the ones specific to the last cell
                // sometimes I may end up with two cell having lines, hence I am doing this twice. There must be a better way,
                // maybe to throttle the mouse event?
                remove_lines()
                // remove_line(last_visited) // I am throttling now, hence this is now redundant
                // $('html,body').css('cursor', 'pointer');
                cellMouseHover(instanceId.label)
                last_visited = instanceId.label
            }
        }
        else {
            $('html,body').css('cursor', 'default');
            clearScreen()
        }
    }
    else {
        // if you are now hovering over any cell, remove any lines you have drawn already
        // and map the last_visited variable to 0 (ie the label for the background)
        remove_lines()
        // $('html,body').css('cursor', 'default');
        last_visited = 0
    }
}


function cellMouseHover(label) {
    console.log('Hovering over cell: ' + label)
    // "https://storage.googleapis.com/merfish_data/cellData/"
    // d3.json("./py/cellData/" + label + ".json", outer(label));
    d3.queue()
        // .defer(d3.json, "https://storage.googleapis.com/izzie_sfn/cellData/"+ label + ".json")
        .defer(d3.json, "data/cell_gene_counts/"+ label + ".json")
        // .defer(d3.json, "./py/cellData/" + label + ".json")
        //.defer(d3.json, "viewer//libs/glyphConfig.")
        .await(splitArgs(label))
}

function splitArgs(label) {
    return (err, ...args) => {

        var data = args[0];
        var targetCell = cellData.filter(d => d.label === label)[0]
        var lines = make_line(data, targetCell)
        lines.map(d => viewer.scene.scene.add(d));
        var spots = groupBy(data, 'gene');
        showControls()
        renderDataTable(spots, targetCell)
        donutchart(targetCell)

    }
}

// function outer(label){
//     return function onCellMouseHover(data) {
//         var targetCell = cellData.filter(d => d.label === label)[0]
//         var lines = make_line(data, targetCell)
//         lines.map(d => viewer.scene.scene.add(d));
//         var spots = groupBy(data, 'gene');
//         // $('#dataTableControl').show();
//         // $('#cellCoordsControl').show();
//         renderDataTable(spots, targetCell)
//         donutchart(targetCell)
//     }
// }

function get_color(gene){
    let hex;
    hex = glyphSettings().filter(d => d.gene === gene)
    if (hex.length === 0) {
       hex = glyphSettings().filter(d => d.gene === 'generic')
    }
    return d3.color(hex[0].color).rgb()
}

function make_line(obj, targetCell, geneColors){
    var arr = Object.entries(obj).map(d => d[1]).flat()
    arr.forEach(d => d['r'] = get_color(d.gene).r)
    arr.forEach(d => d['g'] = get_color(d.gene).g)
    arr.forEach(d => d['b'] = get_color(d.gene).b)
    var out = arr.map(d => {
        return make_line_helper(d, targetCell)
    });
    return out
}

function remove_line(label){
    var scene = viewer.scene.scene
    scene.children.filter(d => (d.type === "Line") && (d.name === label)).forEach(el => scene.remove(el))
}

function make_line_helper(spotData, targetCell) {
    var points = [];
    points.push(
        new THREE.Vector3(spotData.x, spotData.y, spotData.z),
        new THREE.Vector3(targetCell.x, targetCell.y, targetCell.z)
    )
    var geometry = new THREE.BufferGeometry().setFromPoints(points);
    // CREATE THE LINE
    var line = new THREE.Line(
        geometry,
        new THREE.LineBasicMaterial({
            color:  new THREE.Color( spotData.r/255.0, spotData.g/255.0, spotData.b/255.0)
        }),
    );
    line.name = targetCell.label
    return line
}

export default initScene