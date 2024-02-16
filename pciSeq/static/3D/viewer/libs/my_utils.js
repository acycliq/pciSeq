import * as THREE from "./three.js/build/three.module.js";

export function getAllPointsOfPointCloud(pointCloud) {
    var list = [];
    var array = pointCloud.pcoGeometry.root.geometry.attributes.position.array;
    var index = 0;
    for (var i = 0; i < pointCloud.pcoGeometry.root.geometry.attributes.position.count;i=i+3) {
        var x = array[i + 0];
        var y = array[i+ 1];
        var z = array[i + 2];
        let position = new THREE.Vector3(x, y, z);
        position.applyMatrix4(pointCloud.matrixWorld);
        list[index] = position;
        index++;
    }
    return list;
}

export function instanceShow( im, instanceIdx, visible )
{
    // from https://discourse.threejs.org/t/how-to-show-and-hide-an-instance-in-instance-mesh/28198/15
    // im - InstancedMesh
    // instanceIdx - index of an instance
    // visible - boolean show or hide
    const numVisibleInst = im.count;
    const maxInstances = im.instanceMatrix.count;
    var mtx = new THREE.Matrix4();
    var mtx2 = new THREE.Matrix4();

    if( instanceIdx >= maxInstances )
        throw Error("!");

    let lastInstIdx;

    // if show instance which is hidden
    if( visible && instanceIdx >= im.count )
    {
        lastInstIdx = im.count;	// last instance idx
        im.count += 1; 		// show last
    }

    // if hide instance which is visible
    if( !visible && instanceIdx < im.count )
    {
        im.count -= 1; 			// hide last instance
        lastInstIdx = im.count;	// last instance
    }

    // Swap instances: instanceIdx <-> lastInstIdx

    // Swap matrices
    im.getMatrixAt( lastInstIdx, mtx );
    im.getMatrixAt( instanceIdx, mtx2 );

    im.setMatrixAt( lastInstIdx, mtx2 );
    im.setMatrixAt( instanceIdx, mtx );

    im.instanceMatrix.needsUpdate = true;	// !
}

export function hideAll(){
    var count = cells.front_face.instancedMesh.instanceMatrix.count
    for (var i=0; i<count; i++){
        instanceShow(cells.front_face.instancedMesh, i, false)
    }
}

export function showAll(){
    var count = cells.front_face.instancedMesh.instanceMatrix.count
    for (var i=0; i<count; i++){
        instanceShow(cells.front_face.instancedMesh, i, true)
    }
}



