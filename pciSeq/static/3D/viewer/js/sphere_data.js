function get_cell_xyz() {
    var data = [],
        sphere_position = [
            {x: 35, y: 15, z: 17.5},
            {x: -44, y: 23, z: -2.5},
            {x: 43, y: 100, z: -3.5},
            {x: 43, y: 88, z: -15.5},
            {x: 43, y: 126.5, z: 15},
            {x: -5, y: 115.5, z: 7},
            {x: -73, y: -121.5, z: 5.5},
            {x: -66, y: -150.5, z: 17.5},
            {x: -53, y: -27, z: -7.5},
            {x: 10, y: -95, z: -3},
            {x: -20, y: -90, z: -12},
            {x: -13, y: -175, z: -3},
            {x: -125, y: -15, z: 20},
        ],
        sphere_scale = [
            {x: 7.5, y: 10, z: 3.0},
            {x: 15, y: 25, z: 10},
            {x: 13, y: 8, z: 4},
            {x: 12, y: 8, z: 2},
            {x: 8, y: 8, z: 8},
            {x: 12, y: 12, z: 8},
            {x: 17, y: 12, z: 6},
            {x: 12, y: 10, z: 8},
            {x: 14, y: 13, z: 10},
            {x: 13, y: 8, z: 8},
            {x: 8, y: 8, z: 8},
            {x: 12, y: 8, z: 8},
            {x: 8, y: 35, z: 8},
        ],
        sphere_rotation = [
            {x: 0.0, y: 0.0, z: 0.0},
            {x: -Math.PI / 8, y: 0, z: -Math.PI / 4},
            {x: 0, y: 0, z: 0},
            {x: 0, y: 0, z: 0},
            {x: 0, y: 0, z: 0},
            {x: -Math.PI / 8, y: Math.PI / 8, z: 0},
            {x: 0, y: Math.PI / 8, z: Math.PI / 6},
            {x: 0, y: 0, z: 0},
            {x: 0, y: 0, z: 0},
            {x: 0, y: 0, z: 0},
            {x: 0, y: 0, z: 0},
            {x: -Math.PI / 8, y: 0, z: Math.PI / 4},
            {x: 0, y: -Math.PI, z: Math.PI / 10},
        ],
        sphere_color = [
            {r: 0, g: 0, b: 0.9},
            {r: 0, g: 0, b: 0.9},
            {r: 0, g: 0.701960, b: 1},
            {r: 1, g: 0, b: 230 / 255},
            {r: 1, g: 0, b: 230 / 255},
            {r: 0, g: 0.701960, b: 1},
            {r: 0, g: 0.701960, b: 1},
            {r: 1, g: 0, b: 230 / 255},
            {r: 0, g: 0.701960, b: 1},
            {r: 0, g: 1, b: 0},
            {r: 100 / 255, g: 0, b: 90 / 255},
            {r: 100 / 255, g: 0, b: 90 / 255},
            {r: 1, g: 1, b: 1},
        ];

    for (var i = 0; i < sphere_position.length; i++) {
        var dp = sphere_position[i],
            ds = sphere_scale[i],
            dr = sphere_rotation[i],
            dc = sphere_color[i];
        data.push({
            position: new THREE.Vector3(dp.x, dp.y + 16.0, dp.z),
            scale: new THREE.Vector3(ds.x, ds.y, ds.z),
            rotation: new THREE.Vector3(dr.x, dr.y, dr.z),
            color: dc,
        });
    }

    return data
}
