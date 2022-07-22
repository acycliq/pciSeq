function get_sim_cell_xyz(n) {
    var data = [],
        sphere_position = sim_sphere_position(n),
        sphere_scale = sim_sphere_scale(n),
        sphere_rotation = sim_sphere_rotation(n),
        sphere_color = sim_sphere_color(n);

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

function sim_sphere_position(n){
    var out = [];
    for (var i=0; i<n; ++i){
        var _x = CONFIGSETTINGS.img_width * Math.random() - CONFIGSETTINGS.img_width / 2,
            _y = CONFIGSETTINGS.img_height * Math.random() - CONFIGSETTINGS.img_height / 2,
            _z = CONFIGSETTINGS.img_depth * Math.random() - CONFIGSETTINGS.img_depth / 2;
        out.push({x: _x, y: _y, z: _z})
    }
    return out
}

function sim_sphere_scale(n){
    var out = [];
    for (var i=0; i<n; ++i){
        out.push({x: 13 * Math.random(), y: 13* Math.random(), z: 13 * Math.random()})
    }
    return out
}

function sim_sphere_rotation(n){
    var out = [];
    for (var i=0; i<n; ++i){
        out.push({x: Math.random() * Math.PI, y: Math.random() * Math.PI, z: Math.random() * Math.PI})
    }
    return out
}

function sim_sphere_color(n){
    var out = [];
    for (var i=0; i<n; ++i){
        out.push({r: Math.random(), g: Math.random(), b: Math.random()})
    }
    return out
}

