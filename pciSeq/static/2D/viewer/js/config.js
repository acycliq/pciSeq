function config() {
    return {
            roi: {"x0": 0, "x1": 7602, "y0": 0, "y1": 5471 },
            zoomLevels: 10, // maximum zoom levels. Leave that at 10.
            tiles: 'https://storage.googleapis.com/ca1-data/img/262144px/{z}/{y}/{x}.jpg',
            cellData: [{mediaLink: './data/cellData.tsv', size: "2396160"}],
            geneData: [{mediaLink: './data/geneData.tsv', size: "9560621"}],
            cellBoundaries: [{mediaLink: './data/cellBoundaries.tsv', size: "1306209"}],
        }
}
