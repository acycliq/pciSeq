function renderDataTable(d) {

        
    var mydata = [];
    var mydata2 = [];

    // var str = "<strong>Cell Num: </strong>" + d.Cell_Num +
    //     ",  (<strong>x, y</strong>): (" + d.x.toFixed(2) + ", " + d.y.toFixed(2) + ")";
    // document.getElementById('dtTitle').innerHTML = str;
    var n = d3.max([d.CellGeneCount.length, d.Genenames.length]);
    for (i = 0; i < n; i++) {
        mydata.push({
            "Genenames": (d.Genenames[i] === undefined) ? "" : d.Genenames[i],
            "CellGeneCount": (d.CellGeneCount[i] === undefined) ? "" : +d.CellGeneCount[i].toFixed(2),
        })
    }

    var n = d3.max([d.ClassName.length, d.Prob.length]);
    for (i = 0; i < n; i++) {
        mydata2.push({
            "ClassName": (d.ClassName[i] === undefined) ? "" : d.ClassName[i],
            "Prob": (d.ClassName[i] === undefined) ? "" : (!Array.isArray(d.Prob)) ? [d.Prob] : d.Prob[i], // d.Prob can be just a float, make sure it is an array
        })
    }


    // check if a there is a reference to a datatable.
    // If yes, refresh with the new data
    // Otherwise create and populate a datatable
    if ($.fn.dataTable.isDataTable('#dtTable')) {
        table = $('#dtTable').DataTable();
        table.clear().rows.add(mydata).draw();
    } else {
        table = $('#dtTable').DataTable({
            //bFilter: false,
            "lengthChange": false,
            searching: false,
            //"scrollY":        "200px",
            //"scrollCollapse": true,
            "paging": true,
            "pagingType": "simple",
            // "dom": 't',
            "bInfo": false, //Dont display info e.g. "Showing 1 to 4 of 4 entries"
            // "paging": false,//Dont want paging
            "bPaginate": false,//Dont want paging

            "data": mydata,
            "columns": [
                    {
                        title: "Gene",
                        data: "Genenames"
                    },
                    {
                        title: "Counts",
                        data: "CellGeneCount"
                    },
                ],
        });

    }

    function getTotal(table){
        var total = table.column(1).data().reduce(function (a,b) {return a+b}, 0)
        $(table.column(1).footer()).html('Total: ' +total )
        console.log('Total number of gene counts: ' + total )

        return total
    }


    if ($.fn.dataTable.isDataTable('#dtTable2')) {
        table2 = $('#dtTable2').DataTable();
        table2.clear().rows.add(mydata2).draw();
    } else {
        table2 = $('#dtTable2').DataTable({
            //bFilter: false,
            "lengthChange": false,
            searching: false,
            //"scrollY":        "200px",
            //"scrollCollapse": true,
            "paging": true,
            //dom: 't',
            "data": mydata2,
            "columns": [
                {
                    title: "Class Name",
                    data: "ClassName"
                            },
                {
                    title: "Prob",
                    data: "Prob"
                            },
                          ]
        });
    }

    // Sort by column 1 and then re-draw
    table
        .order([1, 'desc'])
        .draw();

    table2
        .order([1, 'desc'])
        .draw();

    var total = getTotal(table);
        // centroid = dapiConfig.t.untransform(d.centroid);

    var str = "<b> <strong>Cell Num: </strong>" + d.cell_id
        + ", <strong>Gene Counts: </strong>" + total.toFixed(0)
        + ",  (<strong>x, y</strong>): (" + d.centroid[0].toFixed(0) + ", " + d.centroid[1].toFixed(0) + ") </b>";

    if (pinnedControls){
        str = str + "<img src='https://cdn2.iconfinder.com/data/icons/oxygen/48x48/actions/note2.png' class='ribbon'/>";
    }
    else
        str = str +  "<img src='https://cdn2.iconfinder.com/data/icons/snipicons/500/pin-128.png' class='ribbon'/>";
    document.getElementById('dtTitle').innerHTML = str

}
