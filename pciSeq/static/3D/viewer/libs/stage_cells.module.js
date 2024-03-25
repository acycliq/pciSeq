import * as THREE from "./three.js/build/three.module.js";

function make_cells_2(data) {
    var front_props = {
            side: THREE.FrontSide,
            opacity: 0.4,
            name: 'front_mesh'
        },
        back_props = {
            side: THREE.BackSide,
            opacity: 0.9,
            name: 'back_mesh'
        };

    // remove the zero class cells. No need to plot them
    var NON_ZERO_CELLS = data.filter(d => d.topClass !== 'ZeroXXX');
    // NON_ZERO_CELLS = NON_ZERO_CELLS.filter(d => d.topClass.startsWith('Calb2'))
    // data = [data[0]]
    var front_face = ellipsoids_2(NON_ZERO_CELLS, front_props),
        back_face = ellipsoids_2(NON_ZERO_CELLS, back_props);
        var cells = {};
    cells.front_face = front_face;
    cells.back_face = back_face;
    return cells
}

function ellipsoids_2(data, props) {
    var counts = data.length,
        loader = new THREE.TextureLoader();

    const flakesTexture = loader.load('./viewer/flakes.png')
    const base_props = {
        clearcoat: 1.0,
        clearcoatRoughness: 0,
        metalness: 0.065,
        roughness: 0.3,
        normalMap: flakesTexture,
        normalScale: new THREE.Vector2(0.3, 0.3),
        transmission: 0.0,
        transparent: true,
        // envMap: true,
    };
    var material = new THREE.MeshPhysicalMaterial(base_props);

    material.side = props.side;
    material.color = props.color;
    material.opacity = props.opacity;
    material.normalMap.wrapS = material.normalMap.wrapT = THREE.RepeatWrapping;
    material.normalMap.repeat = new THREE.Vector2(30, 30)


    var uScale = 0;
    var widthSegments = 16,
        heightSegments = 0.5 * widthSegments;
    var geometry =  new THREE.SphereBufferGeometry(1, widthSegments, heightSegments);
    var _n = geometry.index.count/3;
    console.log('triangles: ' + (_n * counts).toLocaleString());
    var INSTANCEDMESH = new THREE.InstancedMesh(
        //provide geometry
        geometry,

        //provide material
        material,

        //how many instances to allocate
        counts,

        //is the scale known to be uniform, will do less shader work, improperly applying this will result in wrong shading
        !!uScale
    );

    // var dummy = new THREE.Object3D();
    var bbox_items = [];
    console.log('tic')
    var temp_obj = new THREE.Mesh(new THREE.SphereGeometry(), new THREE.MeshBasicMaterial());
    for (var i = 0; i < counts; i++) {
        var coords = data[i],
            scales = {x:data[i].sphere_scale[0], y:data[i].sphere_scale[1], z:data[i].sphere_scale[2]},
            rot = {x:data[i].sphere_rotation[0], y:data[i].sphere_rotation[1], z:data[i].sphere_rotation[2]},
            rgb = {r: data[i].r, g: data[i].g, b: data[i].b};
            // topClass = data[i].topClass;
            // color =  data[i].color;
        var dummy = new THREE.Object3D();
        dummy.position.set(coords.x, coords.y, coords.z);
        dummy.scale.set(scales.x*0.99, scales.y*0.99, scales.z*0.99);
        dummy.rotation.set(rot.x, rot.y, rot.z);
        dummy.updateMatrix();
        INSTANCEDMESH.name = props.name;
        INSTANCEDMESH.setMatrixAt(i, dummy.matrix);
        INSTANCEDMESH.setColorAt(i, new THREE.Color( rgb.r/255.0, rgb.g/255.0, rgb.b/255.0 ));
        temp_obj.applyMatrix4(dummy.matrix)
    }
    console.log('toc')
    INSTANCEDMESH.instanceColor.needsUpdate = true;
    INSTANCEDMESH.visible = true;
    INSTANCEDMESH.castShadow = true;
    INSTANCEDMESH.receiveShadow = false;

    function LOD_ramp() {
        // camera.position.distanceTo(scene.position) < 300? mesh_LOD(): null
        if (camera.position.distanceTo(scene.position) < 300) {
            console.log('Less than 300')
        }

    }


    // viewer.scene.scene.add(new THREE.AmbientLight(0x666666));

    var ellipsoidData = {};
    ellipsoidData.instancedMesh = INSTANCEDMESH;
    ellipsoidData.LOD_ramp = LOD_ramp;


    return ellipsoidData
}

function count_triangles(m){
    // input m is the mesh
    var _n = m.geometry.index.count/3;
    var count = m.count;
    console.log('triangles: ' + (_n * count).toLocaleString());
}


export function tree(data) {
    // makes the tree object to pass into the tree control as an overlay
    var mapper = {},
        root = {
            text: 'Cell Classes',
            selectAllCheckbox: 'Un/select all',
            children: []
        };

    for (var str of data) {
        var sep = '.', //configSettings.class_name_separator,
            splits,
            text = '';
        // let splits = str.match(/[a-zA-Z]+|[0-9]+/g), //str.split('.'),
        if (sep === '') {
            console.log('Assuming that class name is a string followed by a number, like Astro1, Astro2 etc');
            splits = str.match(/[a-zA-Z]+|[0-9]+/g) //str.split('.'),
        } else {
            splits = str.split(sep)
        }
        ;
        splits.reduce(myReducer(text), root)
    }

    function myReducer(text) {
        return function (parent, place, i, arr) {
            if (text) {
                var sep = '.'; //configSettings.class_name_separator;
                text += sep + `${place}`; // `.${place}`;
            } else
                text = place;

            if (!mapper[text]) {
                var o = {text: text};
                // o.id = text;
                o.collapsed = true;
                if (i === arr.length - 1) {
                    o.layer = text
                    o.id = text
                    // o.layer = masterCellContainer.getChildByName(label);
                }
                mapper[text] = o;
                parent.selectAllCheckbox = true;
                parent.children = parent.children || [];
                parent.children.push(o)
            }
            return mapper[text];
        }
    }

    return root
}

export function myjsTree(treeData) {

    // Create an jstree instance
    $('#jstree_demo').jstree({ // config object start

      "core": {                    // core config object
        "mulitple": false,         // disallow multiple selection
        "animation": 100,          // 200ms is default value
        "check_callback" : true,   // this make contextmenu plugin to work
        "themes": {
            "variant": "medium",
            "dots": false,
            "icons":false
        },

        "data": treeData
        //         [
        //   // The required JSON format for populating a tree
        //   { "text": "Root",
        //     "state": {
        //       "opened": true
        //     },
        //    "type": "demo",
        //    "children": [ // "children" key is an array of objects
        //       { "text": "Child node 1",
        //         "li_attr": {
        //           "class": "li-style"
        //        },
        //        "a_attr": {
        //          "class": "a-style"
        //        }
        //       },
        //       { "text": "Child node 2",
        //        "state": {
        //          "opened": true
        //        },
        //        "children": [
        //          { "text": "Grandchild 1",
        //            "state": {
        //              // "disabled": true,
        //              "selected": true
        //           }
        //          },
        //          { "text": "Grandchild 2",
        //            "state": {
        //              "disabled": true
        //           },
        //           "children": [ "Grandson 1", "Grandson 2" ]
        //          },
        //          { "text": "Grandchild 3",
        //            "icon": "glyphicon glyphicon-upload"
        //          },
        //        ]
        //       },
        //     {"text": "Child node 3",
        //      // if all we need is node with the given text, it can be written like this
        //      "children": [ "Grandchild 1", "Grandchild 1" ]
        //     },
        //    ]} // root node end, end of JSON
        // ] // data core options end

      }, // core end

      // Types plugin
      "types" : {
        "default" : {
          "icon" : "glyphicon glyphicon-flash"
        },
        "demo" : {
          "icon" : "glyphicon glyphicon-th-large"
        }
      },

      // config object for Checkbox plugin (declared below at plugins options)
      "checkbox": {
        "keep_selected_style": false,  // default: false
        "three_state": true,           // default: true
        "whole_node": true             // default: true
      },

      "conditionalselect" : function (node, event) {
        return false;
      },

      // injecting plugins
      "plugins" : [
            "checkbox",
            "contextmenu",
            // "dnd",
            // "massload",
            // "search",
            // "sort",
            // "state",
            "types",
            // Unique plugin has no options, it just prevents renaming and moving nodes
            // to a parent, which already contains a node with the same name.
            "unique",
            // "wholerow",
            // "conditionalselect",
            "changed"
      ]
    }); // config object end

     // AJAX loading JSON Example:
     $('#jstree_ajax_demo').jstree({
      'core': {
        'data': {
          "url" : "https://codepen.io/stefanradivojevic/pen/dWLZOb.js",
          "dataType" : "json" // needed only if you do not supply JSON headers
        }
      },
      // Types plugin
      "types" : {
        "default" : {
          "icon" : "glyphicon glyphicon-record"
        }
      },
       "plugins" : [ "types", "unique" ]
    });

    // Listen for events - example
    $('#jstree_demo').on("changed.jstree", function (e, data) {
      // changed.jstree is a event
      console.log('selected: ' + data.changed.selected);
      console.log('deselected: ' + data.changed.deselected);
    });

  };

export default make_cells_2

