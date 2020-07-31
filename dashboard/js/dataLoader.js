function data_loader(workPackage) {
    var previous_avg = 0;
    const average = list => list.reduce((prev, curr) => prev + curr) / list.length;
    const getFileName = (d) => d.split('/').pop();
    const getRootName = (d) => d.split('/').pop().split('_')[0];
    const agg_counts = (d) => {
        var out = {};
        var key;
        d.forEach(function (i) {
            key = getRootName(i);
            out[key] = (out[key] || 0) + 1;
        });
        return out
    };

    const sort_dict = (d) => {
        out = {};
        var keys = Object.keys(d);
        keys.sort();

        for (var i = 0; i < keys.length; i++) {
            out[keys[i]] = d[keys[i]];
        }
        return out
    };

    function setupDict(arr_in) {
        this.counts = agg_counts(arr_in); // returns a dict with keys the names of the datasets and values how many split files each one has
        var key_names = Object.keys(this.counts);  // the dataset name
        this.data = {}; //keeps the data from the tsv flatfiles
        this.perc = {}; //keeps how much of the data have been parsed (in percentage)
        this.bytesStreamed = {}; // Keeps how many bytes in total have been streamed
        this.filePartSize = {};  // keeps the size in bytes of the (split) flatfiles
        this.filePartName = {};  // keeps the (unique) names of the (split) flatfiles already visited
        this.fileSize = {};
        for (var i = 0; i < key_names.length; i++) {
            this.data[key_names[i]] = [];
            this.perc[key_names[i]] = 0;
            this.bytesStreamed[key_names[i]] = 0;
            this.filePartSize[key_names[i]] = [];
            this.filePartName[key_names[i]] = [];
        }

        this.update = (event) => {
            if (event.data && event.data.url) {
                var key = getRootName(event.data.url);
                this.data[key] = this.data[key].concat(event.data.items);
                this.bytesStreamed[key] += event.data.byteStats[3];
                var filename = getFileName(event.data.url);
                if ( !this.filePartName[key].includes(filename) ) {
                    // if you stream the file for the very first time, add it to the filePartName dict
                    // and update the filePartSize dict with its size.
                    this.filePartName[key] = this.filePartName[key].concat([filename]);
                    var fSize = +event.data.byteStats[2];
                    this.filePartSize[key] = this.filePartSize[key].concat([fSize]);
                }
                var totalSize = +this.get_filesize(key);
                var x = +this.bytesStreamed[key] / totalSize;
                this.perc[key] = Math.floor(x * 100);
            }
        }

        this.get_filesize = (key) => {
            // if you have visited all the split files for the given key then calc the total size.
            // Otherwise give an estimate.
            var N = this.counts[key],
                arr = this.filePartSize[key];
            if (arr.length != N){
                return arr[0] * N
            }
            else {
                // sum all the elements of the array
                return arr.reduce((a, b) => a + b, 0)
            }

        }
    }


    function makeInnerHtml(data) {
        var innerHtml = "";
        for (var i = 0; i < data.length; i++) {
            innerHtml = innerHtml + data[i] + "/"
        }
        innerHtml = innerHtml.slice(0, -1); //remove the last slash from the end of the string
        innerHtml = "<span> Loaded: " + innerHtml + "</span>";

        return innerHtml
    }


    function setupWorker() {
        // create a web worker that streams the chart data
        worker = new Worker("./streaming-tsv-parser.js");
        worker.onmessage = function (event) {
            stats.update(event);
            if (event.data.finished) {
                // console.log(data)
                document.getElementById("done").innerHTML = "<span>All Done! </span>";
                onDataLoaded(stats.data);
                // redraw(stats);
                return
            }
            redraw(stats);
        };
    }

    setupWorker();
    var stats = new setupDict(workPackage);
    worker.postMessage(workPackage);

    function redraw(stats) {

        var data = stats.data;
        var perc = stats.perc;

        innerHtml = makeInnerHtml(Object.values(data).map(d => d.length));
        document.getElementById("loading").innerHTML = innerHtml;

        innerHtml = makeInnerHtml(Object.values(perc).map(d => d + '% '));
        innerHtml = innerHtml.replace(/% \//g, '% ');
        document.getElementById("loading_perc").innerHTML = innerHtml;

        var perc_array = Object.values(perc);
        var min_perc = Math.min(...perc_array);

        if (average(perc_array) >= previous_avg + 2){
            // refresh the progress bar animation, incrementally, every 2%
            bar.animate(average(perc_array) / 100); // Number from 0.0 to 1.0
            previous_avg = average(perc_array);
        }

    }

    function onDataLoaded(data){
        [cellBoundaries, cellData] = postLoad([data.cellCoords, data.cellData]);

        // sort cellBoundaries and cellData. These two should be aligned
        cellBoundaries.sort((a, b) => a.cell_id - b.cell_id);
        cellData.sort((a, b) => a.cell_id - b.cell_id);

        all_geneData = data.geneData;
        console.log('loading data finished');
        console.log('num of genes loaded: ' + all_geneData.length);
        console.log('num of cells loaded: ' + cellData.length);


        //finaly make a spatial index on the spots. We will need that to filter them if/when needed
        console.log('Creating the index');
        spotsIndex = new KDBush(all_geneData, p => p.x, p => p.y, 64, Int32Array);

        // do now the chart
        dapiChart(configSettings);

    }
}
