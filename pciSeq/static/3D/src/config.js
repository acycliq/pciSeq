function config() {
    var ini = [
        {   // 1.
            name: 'default',
            img_width: 5905,
            img_height: 5882,
            img_depth: 15 * 6.0121,
            particle_size: 8000.0,
            // zThres: 1500.0,
            geneData: 'https://api.github.com/repos/acycliq/B2A3/contents/data/geneData?ref=master',
            cellData: 'https://api.github.com/repos/acycliq/B2A3/contents/data/cellData?ref=master',
        },
    ];
    return d3.map(ini, function (d) {return d.name;})
}

