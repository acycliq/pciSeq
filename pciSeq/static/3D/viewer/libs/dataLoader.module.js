import make_cells_2 from './stage_cells.module.js'
import initScene from './initScene.module.js'

function data_loader(workPackage) {
    var data = [],
        agg_data,
        previous_avg = 0;
    var average = list => list.reduce((prev, curr) => prev + curr) / list.length;
    workPackage = workPackage.sort((a,b) =>  a.size-b.size); //order by size (hemmm...doest really matter, does it?? Everything happens in parallel)

    function aggregate_stats(workPackage){
        var out = [];
        var names = [... new Set(workPackage.map(d => d.root_name))].sort();
        names.forEach(name => {
            var temp = workPackage.filter(d => d.root_name === name);
            // console.log(temp)
            var size = temp.map(d => d.size).reduce((a, b) => a + b, 0);
            var bytes_streamed =  temp.map(d => d.bytes_streamed).reduce((a, b) => a + b, 0);
            var data_length =  temp.map(d => d.data_length).reduce((a, b) => a + b, 0);
            var progress = bytes_streamed / size;
            out.push({name, size, bytes_streamed, progress, data_length})
            // console.log(out)
        });
        return out
    }

    function aggregate_data(workPackage) {
        var out = {};
        var names = [...new Set(workPackage.map(d => d.root_name))].sort();
        names.forEach(name => {
            var temp = workPackage.filter(d => d.root_name === name);
            var data = temp.map(d => d.data).reduce((a,b)=>a.concat(b));
            out[name] = data
        });
        return out
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

    var worker
    function setupWorker() {
        // create a web worker that streams the chart data
        worker = new Worker("./viewer/libs/streaming-tsv-parser.js");
        worker.onmessage = function (event) {
            if (event.data.finished) {
                console.log(agg_data);
                data = aggregate_data(workPackage);
                onDataLoaded(data);
                // redraw(stats);
                return
            }
            var i = event.data.i;
            workPackage[i].bytes_streamed += event.data.bytes_streamed;
            workPackage[i].data = workPackage[i].data.concat(event.data.items);
            workPackage[i].data_length += event.data.items.length;
            agg_data = aggregate_stats(workPackage);

            redraw(agg_data);
        };
    }

    function redraw(data) {
        // console.log(data[0])
        // console.log(data[1])
        var innerHtml = makeInnerHtml(data.map(d => d3.format(",")(d.data_length)));
        document.getElementById("loading").innerHTML = innerHtml;

        innerHtml = makeInnerHtml(data.map(d => d3.format(".0%")(d.progress) ));
        innerHtml = innerHtml.replace(/% \//g, '% ');
        document.getElementById("loading_perc").innerHTML = innerHtml;

        innerHtml = makeInnerHtml(data.map(d => (d.bytes_streamed/(1024*1024)).toFixed() + 'MB' ));
        document.getElementById("loading_mb").innerHTML = innerHtml;

        var avg = average(data.map(d => d.progress));
        var avg_mb = average(data.map(d => (d.bytes_streamed/(1024*1024)).toFixed() ));
        var progress_1 = data[0].progress;
            // progress_2 = data[1].progress;

        var inc = 0.0; // controls how often it will be update. Every 2% set inc = 0.02
        if (avg >= Math.min(1, previous_avg + inc)) {
            if (avg > 0.99){
                $('#wait_chart').show();
            }
            // refresh the progress animation
            updateDonutChart('#specificChart', progress_1*100, true);
            var mb_1 = (data[0].bytes_streamed/(1024*1024)).toFixed();
            $('#mb').html(mb_1 + 'MB');
            $('#datapoints').html(d3.format(",")(data[0].data_length));

            // updateDonutChart('#specificChart2', progress_2*100, true);
            // $('#mb2').html((data[1].bytes_streamed/(1024*1024)).toFixed() + 'MB');
            // $('#datapoints2').html(d3.format(",")(data[1].data_length));

            previous_avg = avg;
        }
    }

    function onDataLoaded(data) {
        console.log(data)
        cellData = data.cellData
        cellData = postLoad(cellData)
        initScene(cellData)

        console.log('Done!')
    }

    function postLoad(data){
        // loops over the data and attaches the class with the highest prob
        data.forEach(function(d,i) {
            if (d.class_prob.length>0){
                var imax = d.class_prob.indexOf(Math.max(...d.class_prob));
                var top_class = d.ClassName[imax]
                if (top_class !== "undefined") {
                    this[i]['top_class'] = top_class.toString()
                }
                else{
                    this[i]['top_class'] = top_class
                }

            }
        }, data)

        return data
    }

    setupWorker();
    worker.postMessage(workPackage);
}

export default data_loader