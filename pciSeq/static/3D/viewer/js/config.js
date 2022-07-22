function config() {
    var ini = [
        {   // 1.
            name: 'default',
            img_width: 5905,
            img_height: 5882,
            img_depth: 15 * 6.0121,
            particle_size: 8000.0,
            // zThres: 1500.0,
            geneData: [{mediaLink: '../../data/geneData.tsv', size: "65755002"}],
            cellData: [{mediaLink: '../../data/cellData.tsv', size: "10554015"}],
        },
    ];
    return d3.map(ini, function (d) {return d.name;})
}