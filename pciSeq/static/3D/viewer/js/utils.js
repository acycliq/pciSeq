//
// Some useful common functions
//

const groupBy = (array, key) => {
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

function getTaxonomy(gene) {
    return glyphSettings().filter(d => d.gene === gene)[0].taxonomy
}

function getGlyph(gene) {
    console.log(gene)
    return glyphSettings().filter(d => d.gene === gene)[0].glyphName
}

function getColor(gene) {
    const t = getTaxonomy(gene);
    return glyphColor(t)
}

function getShader(glyph) {
    return glyph === 'star6' ? fShader_star6 :
        glyph === 'star5' ? fShader_star5 :
            glyph === 'diamond' ? fShader_diamond :
                glyph === 'square' ? fShader_square :
                    glyph === 'triangleUp' ? fShader_triangleUp :
                        glyph === 'triangleDown' ? fShader_triangleDown :
                            glyph === 'triangleRight' ? fShader_triangleRight :
                                glyph === 'triangleLeft' ? fShader_triangleLeft :
                                    glyph === 'cross' ? fShader_cross :
                                        glyph === 'plus' ? fShader_plus :
                                            glyph === 'asterisk' ? fShader_asterisk :
                                                glyph === 'circle' ? fShader_circle :
                                                    glyph === 'point' ? fShader_point :
                                                        console.log('Shader for ' + glyph + ' is missing.')
}

function hexToRgb(hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}


function getGenePanel(geneData) {
    var panel = [...new Set(geneData.map(d => d.Gene))].sort();

    // drop a warning if a gene is not set in the configuration file
    var cfg_genes = glyphSettings().map(d => d.gene).sort();
    var missing = panel.filter(x => !cfg_genes.includes(x));
    if (missing.length > 0) {
        console.log('Waring: These genes have not been assigned color, glyph etc in the glyphConfig.js: ' + missing);
    }

    // save the gene panel to the local storage
    sessionStorage.setItem('gene_panel', JSON.stringify(panel));
    console.log('Gene panel written to local storage')

    return panel
}

function legendControl() {
    var legendLink = document.querySelector(`#legend-link`);

    if (!legend_added) {
        legendLink.addEventListener(`click`, () => {
            // Opens the page and stores the opened window instance
            legendWindow = window.open(`./src/genes_datatable.html`); // <--- THAT NEEDS TO BE EXPOSED TO THE USER. MOVE I INSIDE config.js MAYBE
        });
    }
    legend_added = true;

    $('#legend').show()
    console.log('legend added')
}

function hoverPieces() {
    raycaster.setFromCamera(mouse, camera);
    raycaster.far = configSettings.zThres;
    raycaster.params.Points.threshold = 0.5;
    const intersects = raycaster.intersectObjects(scene.children);
    if (intersects.length) {
        console.log('Intersects')
        if (intersected != intersects[0].object) {
            var index = intersects[0].index;
            intersected = intersects[0].object;
            var name = intersected.name,
                uid = name + "_" + index;
            reset_scale(intersected)
            highlight(intersected, name, index);
            console.log('mouseover Gene: ' + name + ' index: ' + index);
            console.log(uid)
            $('html,body').css('cursor', 'pointer');
        }
    } else {
        console.log('does not intersect')
        if (intersected) {
            $('html,body').css('cursor', 'default');
            reset_scale(intersected)
            intersected = null
        }
    }

    // for (let i = 0; i < intersects.length; i++) {
    //     intersects[i].object.material.transparent = true;
    //     intersects[i].object.material.opacity = 0.5;
    // }
}

function highlight(intersected, name, index) {
    if (intersected.type === 'Points') {
        var scales = intersected.geometry.attributes.scale;
        scales.array[index] = 2.0;
        scales.needsUpdate = true;
        console.log('scale ' + name + '_' + index + ' is: ' + intersected.geometry.attributes.scale.array[index])
    }
}

function reset_scale(intersected) {
    if (intersected && intersected.type === 'Points') {
        var scales = intersected.geometry.attributes.scale;
        console.log('gene is: ' + intersected.name)
        for (var i = 0; i < scales.count; i++) {
            if (scales.array[i] != 1.0) {
                console.log('resetting at position: ' + i)
                scales.array[i] = 1.0;
            }
        }
        scales.needsUpdate = true;
    }
}


// window.addEventListener('mousemove', (e) => {
//     var vec = new THREE.Vector3(); // create once and reuse
//     var pos = new THREE.Vector3(); // create once and reuse
//
//     vec.set(
//         ( event.clientX / window.innerWidth ) * 2 - 1,
//         - ( event.clientY / window.innerHeight ) * 2 + 1,
//         0.5 );
//
//     vec.unproject( camera );
//
//     vec.sub( camera.position ).normalize();
//
//     var distance = - camera.position.z / vec.z;
//
//     pos.copy( camera.position ).add( vec.multiplyScalar( distance ) );
//     console.log('mouse x: ' + pos.x)
//     console.log('mouse y: ' + pos.y)
//
// });


function hideAxes() {
    for (var i = 0; i < scene.children.length; i++) {
        if (scene.children[i].name === 'xyz_axes') {
            scene.children[i].visible = false;
            console.log('Axes are switched off.')
        }
    }
}


function showAxes() {
    for (var i = 0; i < scene.children.length; i++) {
        if (scene.children[i].name === 'xyz_axes') {
            scene.children[i].visible = true;
            console.log('Axes are switched on.')
        }
    }
}


function createAxes(axisLength, target) {
    function v(x, y, z) {
        return new THREE.Vector3(x, y, z);
    }

    function axis(p1, p2, color, target) {
        var line;
        var points = [];
        points.push(p1);
        points.push(p2);
        var lineGeometry = new THREE.BufferGeometry().setFromPoints(points);
        // lineGeometry = new THREE.Geometry(),
        lineMat = new THREE.LineBasicMaterial({color: color, lineWidth: 1});
        // lineGeometry.vertices.push( p1, p2 );
        line = new THREE.Line(lineGeometry, lineMat);
        return line
        // target.add( line );
    }

    var group = new THREE.Group();
    group.name = 'xyz_axes'
    var axis_x = axis(v(-axisLength, 0, 0), v(axisLength, 0, 0), 0xFF0000, target);
    var axis_y = axis(v(0, -axisLength, 0), v(0, axisLength, 0), 0x00FF00, target);
    var axis_z = axis(v(0, 0, -axisLength), v(0, 0, axisLength), 0x0000FF, target);
    group.add(axis_x);
    group.add(axis_y);
    group.add(axis_z);

    return group
}

function LOD_ramp(x) {
    // inpout x is camera.position.distanceTo(scene.position)
    var lod = {};
    if (x > 1000) {
        lod.level = 'level_0';
        lod.w = 4;
        lod.h = 4;
    } else if (x > 500) {
        lod.level = 'level_1';
        lod.w = 6;
        lod.h = 4;
    } else if (x > 300) {
        lod.level = 'level_2';
        lod.w = 12;
        lod.h = 6;
    } else if (x > 100) {
        lod.level = 'level_3';
        lod.w = 18;
        lod.h = 8;
        if (SCENE.children.filter(d => d.name === 'back_mesh').length) {
            SCENE.children.filter(d => d.name === 'back_mesh').forEach(d => {
                SCENE.remove(d);
                console.log('back_mesh removed')
            })
        }
    } else {
        lod.level = 'level_4';
        lod.w = 24;
        lod.h = 12;
        if (!SCENE.children.filter(d => d.name === 'back_mesh').length) {
            SCENE.add(INSTANCEDMESH.back_face.instancedMesh);
            console.info('back face added')
        }
    }
    return lod
}

function mesh_LOD(w, h) {
    var meshes = SCENE.children.filter(v => v.type === 'Mesh');
    meshes.forEach(d => {
        const clonedGeometry = new THREE.SphereBufferGeometry(1, w, h);
        // clonedGeometry.parameters.widthSegments = paramsGUI.widthSegments;
        d.geometry.dispose();
        d.geometry = clonedGeometry;

        count_triangles(d)
    })
}


function simulate_spots(counts) {
    var nG = GENEPANEL.length;
    var N = Math.ceil(counts / nG);
    var _sim_data = [];
    for (var i = 0; i < GENEPANEL.length; i++) {
        var temp = new Float32Array(N * 3);
        for (var j = 0; j < 3 * N; j++) {
            var cur = 3 * j;
            temp[3 * j] = CONFIGSETTINGS.img_width * Math.random() - CONFIGSETTINGS.img_width / 2;
            temp[3 * j + 1] = CONFIGSETTINGS.img_height * Math.random() - CONFIGSETTINGS.img_height / 2;
            temp[3 * j + 2] = CONFIGSETTINGS.img_depth * Math.random() - CONFIGSETTINGS.img_depth / 2;
        }
        _sim_data[i] = temp
    }
    return _sim_data
}


// from https://github.com/jonathantneal/convert-colors
function rgb2hex(rgbR, rgbG, rgbB) {
	return `#${((1 << 24) + (Math.round(rgbR * 255 / 100) << 16) + (Math.round(rgbG * 255 / 100) << 8) + Math.round(rgbB * 255 / 100)).toString(16).slice(1)}`;
}

function hsl2rgb(hslH, hslS, hslL) {
	// calcuate t2
	const t2 = hslL <= 50 ? hslL * (hslS + 100) / 100 : hslL + hslS - hslL * hslS / 100;

	// calcuate t1
	const t1 = hslL * 2 - t2;

	// calculate rgb
	const [ rgbR, rgbG, rgbB ] = [
		hue2rgb(t1, t2, hslH + 120),
		hue2rgb(t1, t2, hslH),
		hue2rgb(t1, t2, hslH - 120)
	];

	return [ rgbR, rgbG, rgbB ];
}

function hue2rgb(t1, t2, hue) {
	// calculate the ranged hue
	const rhue = hue < 0 ? hue + 360 : hue > 360 ? hue - 360 : hue;

	// calculate the rgb value
	const rgb = rhue * 6 < 360
		? t1 + (t2 - t1) * rhue / 60
	: rhue * 2 < 360
		? t2
	: rhue * 3 < 720
		? t1 + (t2 - t1) * (240 - rhue) / 60
	: t1;

	return rgb;
}

function hsv2hex(hsvH, hsvS, hsvV) {
    // takes in values between 0 and 1
	return _hsv2hex(hsvH*100, hsvS*100, hsvV*100)
}

function _hsv2hex(hsvH, hsvS, hsvV) {
	return rgb2hex(...hsl2rgb(hsvH, hsvS, hsvV));
}

// function hide_back_face() {
//     console.log('Hiding back face')
//     scene.children.forEach(d => {
//         if ((d.type === 'Mesh') & (d.name === 'back_mesh')) {
//             d.visible = false
//         }
//     })
// }

// function scale_helper(points) {
//     var scales = points.geometry.attributes.scale;
//     var count = scales.count;
//
//     for (var i = 0; i < count; i++) {
//         // dynamically change alphas
//         scales.array[i] *= 0.95;
//         if (scales.array[i] < 0.10) {
//             console.log('dynamically change alphas')
//             scales.array[i] = 1.0;
//         }
//     }
//     scales.needsUpdate = true; // important!
// }

function tree(data) {
    // makes the tree object to pass into the tree control as an overlay
    var mapper = {},
        root = {
            text: 'Cell Classes',
            children: [],
            icon: false,
        };

    for (var str of data) {
        var sep = '.',
            splits,
            text = '';
        // let splits = str.match(/[a-zA-Z]+|[0-9]+/g), //str.split('.'),
        if (sep === '') {
            console.log('Assuming that class name is a string followed by a number, like Astro1, Astro2 etc');
            splits = str.match(/[a-zA-Z]+|[0-9]+/g) //str.split('.'),
        } else {
            splits = str.split(sep)
        };
        splits.reduce(myReducer(text), root)
    }

    function myReducer(label) {
        return function (parent, place, i, arr) {
            if (label) {
                var sep = '.';
                label += sep + `${place}`; // `.${place}`;
            } else
                label = place;

            if (!mapper[label]) {
                var o = {text: label};
                o.opened = false;
                o.icon = false;
                if (i === arr.length - 1) {
                    // o.layer = cellContainer_array.filter(d => d.name === label)[0]
                }
                mapper[label] = o;
                parent.opened = false;
                parent.icon = false,
                parent.children = parent.children || [];
                parent.children.push(o)
            }
            return mapper[label];
        }
    }

    return root
}

function class_to_color_helper(x){
    return classColorsCodes().filter(el => el.className.startsWith(x)).map(d => d.color)
}

function class_to_color(x) {
    var out = [...new Set(x.map(d => class_to_color_helper(d))[0])];
    out = [...new Set(out)]
    return out
}


function toggle_background_spots(){
    console.log("pressed")
    if (document.querySelector('#exampleCheck1:checked') !== null) {
        hide_background_spots()
    }
    else{
        show_background_spots()
    }
}

function hide_background_spots() {
    // loop over the genes
    var names = Object.keys(SPOT_COUNTS)
    for (var i = 0; i < names.length; ++i) {
        var name = names[i]
        var counts = SPOT_COUNTS[name]["non_background"]
        SCENE.children.filter(d => d.name === name)[0].geometry.attributes.position.count = counts

    }
}

function show_background_spots(){
    // loop over the genes
    var names = Object.keys(SPOT_COUNTS)
    for (var i=0; i<names.length; ++i){
        var name = names[i]
        var counts = SPOT_COUNTS[name]["total"]
        SCENE.children.filter(d => d.name === name)[0].geometry.attributes.position.count = counts

    }

}