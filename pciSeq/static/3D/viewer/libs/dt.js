function renderDataTable(spots, cell) {

        
    var mydata = [];
    var mydata2 = [];

    for (gene_name in spots){
        mydata.push({
            "Genenames": gene_name,
            "CellGeneCount": +spots[gene_name].length.toFixed(2),
        })
    }

    var n = cell.class_prob.length;
    for (i = 0; i < n; i++) {
        mydata2.push({
            "ClassName": cell.classes[i],
            "Prob": +cell.class_prob[i], // d.Prob can be just a float, make sure it is an array
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
            "lengthChange": false,
            searching: false,
            "paging": true,
            "bInfo": false, //Dont display info e.g. "Showing 1 to 4 of 4 entries"
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
            "lengthChange": false,
            searching: false,
            "paging": true,
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
   //      // centroid = dapiConfig.t.untransform(d.centroid);
   //
    var str = "<b> <strong>Cell Num: </strong>" + cell.label
        + " <br> <strong>Gene Counts: </strong>" + total.toFixed(0)
        // + " <br>  (<strong>x, y, z</strong>): (" + d.X.toFixed(0) + ", " + d.Y.toFixed(0) + ", " + d.Z.toFixed(0) + ")"
        + " <br>  <strong>x:</strong> " + cell.x.toFixed(0)
        + " <br>  <strong>y:</strong> " + cell.y.toFixed(0)
        + " <br>  <strong>z:</strong> " + cell.z.toFixed(0)
        + " </b>";


   document.getElementById('cellCoordsControlText').innerHTML = str

}
