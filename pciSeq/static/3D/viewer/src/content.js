function iniContent(spots, cellData) {
    console.log('Init Viewer');

    var points = GENEPANEL.map((d, i) => my_particles(spots[i], d));
    points.map(d => SCENE.add(d));

    // SCENE.add(my_particles(spots[5], 'Cck'));
    // SCENE.add(my_particles(spots[21], 'Gad1'));
    // SCENE.add(my_particles(spots[0], 'Aldoc'));

    // // add_spheres();
    // INSTANCEDMESH = make_cells(cellData);
    // SCENE.add(INSTANCEDMESH.back_face.instancedMesh); // needs to be rendered first
    // SCENE.add(INSTANCEDMESH.front_face.instancedMesh);

    INSTANCEDMESH = make_cells_2(cellData);
    SCENE.add(INSTANCEDMESH.front_face.instancedMesh);

    CELL_MESH = SCENE.children.filter(d => d.name==='front_mesh')[0]

    // load_dapi()

}
