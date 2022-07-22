function donut(){

	// var width = 280,
	// 	height = 200,
	// 	radius = Math.min(width, height) / 2;

	var cornerRadius = 3, // sets how rounded the corners are on each slice
        padAngle = 0.015; // effectively dictates the gap between slices

    var width = +d3.select("#pie").select("svg").attr("width"),
        height = +d3.select("#pie").select("svg").attr("height")

    var radius = Math.min(width, height) / 2;

    var svg = d3.select("#pie")
        .select("svg")
        .attr("width", width)
		.attr("height", height)
		.append("g")
		.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")"); // Moving the center point

    svg.append("defs").append("pattern")
        .attr('id','myPattern')
        .attr("width", 4)
        .attr("height", 4)
        .attr('patternUnits',"userSpaceOnUse")
        .append('path')
        .attr('fill','none')
        .attr('stroke','#335553')
        .attr('stroke-width','1')
        .attr('d','M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2' );

	svg.append("g")
		.attr("class", "slices");
	svg.append("g")
		.attr("class", "sliceLabels");
	svg.append("g")
		.attr("class", "lines");

	var pie = d3.pie()
		.sort(null)
		.value(function(d) {
			return d.value;
		});

	var arc = d3.arc()
		.outerRadius(radius * 0.8)
		.innerRadius(radius * 0.4)
		.cornerRadius(cornerRadius)
        .padAngle(padAngle);

	var outerArc = d3.arc()
		.innerRadius(radius * 0.9)
		.outerRadius(radius * 0.9);

	//svg.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

	var key = function(d){ return d.data.label; };
	// var colors = ["#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"];
    // var colors = d3.schemeCategory20
    var colorRamp = classColorsCodes()
    var colorMap = d3.map(colorRamp, function(d) { return d.className; });

    var div = d3.select("body").append("div")
        .attr("class", "toolTip")
        .style('opacity', 0);

	var donutData = {};
    donutData.radius = radius;
	donutData.pie = pie;
    donutData.arc = arc;
    donutData.outerArc = outerArc;
    donutData.key = key;
    donutData.colorMap = colorMap;
    donutData.div = div;
    donutData.svg = svg;


    return donutData;

}

// Donut chart for the popup on the cells
function donutPopup(d){
    var feature = d.feature;
    var data = []
    for (var i=0; i < d.feature.properties.ClassName.length; i++) {
        data.push({
            value: d.feature.properties.Prob[i],
            label: d.feature.properties.ClassName[i],
        })
    }

    // Sort now descreasing order. Maybe encapsulate since the same 
    // piece of code is used in scatter.js before the call to barchart and
    // donutchart.
    // For small values assign it to a separate class labeled 'Other'
    var sdata = [];
    var ClassName;
    for (var i = 0; i < d.feature.properties.ClassName.length; i++) {
        d.feature.properties.Prob[i] < 0.02? ClassName = 'Other': ClassName = d.feature.properties.ClassName[i]
        sdata.push({
            Prob: d.feature.properties.Prob[i],
            labels: ClassName,
        })
    }

    // group by class name. This will eventually sum all the values that fall in the same class.
    //Hence if there is a class 'Other' then is will be assigned with the grand total
    sdata =  d3.nest().key(function(d){
        return d.labels; })
        .rollup(function(leaves){
            return d3.sum(leaves, function(d){
                return d.Prob;})
        }).entries(sdata)
        .map(function(d){
            return { label: d.key, value: d.value};
        });

    // sort in decreasing order
    sdata.sort(function(x, y){
        return d3.ascending(y.value, x.value);
    })

    //overwrite data
    data = sdata

    var width = 40;
    var height = 40;
    var radius = dapiConfig.getRadius(0.5 * (Math.min(width, height) * feature.properties.GeneCountTotal));
    if (2*radius >= Math.max(width, height)){
        // readjust the div size if the pie marker is bigger
        var padding = 0;
        width = 2 * radius + padding;
        height = 2 * radius + padding;
    }
    // I think the centre of the piemarkers is a bit off....maybe for a couple of pixels

    var cornerRadius = 3, // sets how rounded the corners are on each slice
        padAngle = 0.015; // effectively dictates the gap between slices

    var div = d3.create("div")
        .attr("class", "piePopupCustom")
    var svg = div.append("svg")
        .attr("width", width)
        .attr("height", height)
        .append("g")
        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")"); // Moving the center point

    svg.append("g")
        .attr("class", "slices");
    svg.append("g")
        .attr("class", "sliceLabels");
    svg.append("g")
        .attr("class", "lines");

    var pie = d3.pie()
        .sort(null)
        .value(function(d) {
            return d.value;
        });

    var arc = d3.arc()
        .outerRadius(radius * 1.0)
        .innerRadius(radius * 0.0)
        .cornerRadius(cornerRadius)
        .padAngle(padAngle);

    var outerArc = d3.arc()
        .innerRadius(radius * 0.9)
        .outerRadius(radius * 0.9);

    var key = function(d){ return d.data.label; };
    // var colors = ["#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"];
    var colors = d3.schemeCategory20;

    if (d3.select("body").select(".popupToolTip").empty()) {
        var map = d3.select('#mymap')
        var divTooltip = L.DomUtil.create('div', 'popupToolTip', map.node());

        var popupToolTipStyle = "position: absolute; " +
            "z-index: 1000; " +
            "display: none; " +
            "width: auto; " +
            "height: auto; " +
            "background: none repeat scroll 0 0 white; " +
            "border: 0 none; " +
            "border-radius: 8px 8px 8px 8px; " +
            "box-shadow: -3px 3px 15px #888888;  " +
            "color: black; " +
            "font: 12px sans-serif; " +
            "padding: 5px; " +
            "text-align: center;"

        divTooltip.setAttribute("style", popupToolTipStyle);

        //var divTooltip = d3.select("body").append("div").attr("class", "popupToolTip");
    }
    divTooltip = d3.select(".popupToolTip");
    // divTooltip.node().setAttribute("style", "position: absolute; z-index: 1000; display: none; width: auto; height: auto; background: none repeat scroll 0 0 white; border: 0 none; border-radius: 8px 8px 8px 8px; box-shadow: -3px 3px 15px #888888;  color: red; font: 12px sans-serif; padding: 5px; text-align: center;")
    // divTooltip.node().setAttribute("style", popupToolTipStyle);
    var percentFormat = d3.format('.2%');


    var labels = d3.map(data, function (d) {return d.label;}).keys();

    var color = d3.scaleOrdinal()
        .domain(labels)
        .range(colors);

    var colorRamp = classColorsCodes()
    var colorMap = d3.map(colorRamp, function(d) { return d.className; });

    /* ------- PIE SLICES -------*/
    var slice = svg.select(".slices").selectAll("path.slice")
        .data(pie(data), key);

    slice.enter()
        .insert("path")
        .attr("class", "slice")
        .on("mousemove", mousemoveHandler)
        .on("mouseout", function (d) {
            divTooltip.style("display", "none");
        })
        .merge(slice)
        .style("fill", function(d) { return colorMap.get(d.data.label).color; })
        .transition().duration(1000)
        .attrTween("d", function(d) {
            this._current = this._current || d;
            var interpolate = d3.interpolate(this._current, d);
            this._current = interpolate(0);
            return function(t) {
                return arc(interpolate(t));
            };
        });

    // svg.append("g")
    //     .append('text')
    //     .attr('class', 'toolCircle')
    //     .attr('dy', 5) // hard-coded. can adjust this to adjust text vertical alignment in tooltip
    //     .text(Math.round( d.feature.properties.GeneCountTotal )) // add text to the circle.
    //     .style('font-size', 0.7*radius+"px")
    //     .style('font', 'sans-serif')
    //     .style('font-weight', 'bold')
    //     .style('text-anchor', 'middle'); // centres text in tooltip


    function mousemoveHandler() {
        divTooltip.style("left", d3.event.pageX + 10 + "px");
        divTooltip.style("top", d3.event.pageY - 25 + "px");
        divTooltip.style("display", "inline-block");
        divTooltip.html((this.__data__.data.label) + "<br>" + percentFormat(this.__data__.data.value) );
    }
        

return div.node();
}

function donutchart(dataset) {
    var percentFormat = d3.format('.2%');

    var data = []
    for (var i=0; i < dataset.ClassName.length; i++) {
        data.push({
            // value: Math.floor(dataset[i].value*10000)/100,
            // label: dataset[i].label,
            value: dataset.Prob[i],
            label: dataset.ClassName[i],
        })
    }


    // For small values assign it to a separate class labeled 'Other'
    // ok, that can be simplified, do this inside the loop above maybe?
    var sdata = [];
    var ClassName;
    for (var i = 0; i < dataset.Prob.length; i++) {
        dataset.Prob[i] < 0.02? ClassName = 'Other': ClassName = dataset.ClassName[i]
        sdata.push({
            Prob: dataset.Prob[i],
            labels: ClassName,
        })
    }

    // group by class name. This will eventually sum all the values that fall in the same class.
    //Hence if there is a class 'Other' then is will be assigned with the grand total
    sdata =  d3.nest().key(function(d){
        return d.labels; })
        .rollup(function(leaves){
            return d3.sum(leaves, function(d){
                return d.Prob;})
        }).entries(sdata)
        .map(function(d){
            return { label: d.key, value: d.value};
        });

    // sort in decreasing order
    sdata.sort(function(x, y){
        return d3.ascending(y.value, x.value);
    })

    //overwrite data
    data = sdata


    var svg = d3.select("#pie").select("svg")
    if (svg.select('.slices').empty()) {
    	donutData = donut();
    }

    var labels = d3.map(data, function (d) {return d.label;}).keys();
    
    // var color = d3.scaleOrdinal()
	// .domain(labels)
	// .range(donutData.colors);
  
	/* ------- PIE SLICES -------*/
	var slice = svg.select(".slices").selectAll("path.slice")
		.data(donutData.pie(data), donutData.key);

	slice.enter()
		.insert("path")
		.attr("class", "slice")
		.on("mousemove", mousemoveHandler)
        .on("mouseout", function (d) {
            donutData.div.style("display", "none");
        })
    .merge(slice)
        //.style("fill", 'url(#myPattern)')
        .style("fill", function(d) { return donutData.colorMap.get(d.data.label).color; })
		.transition().duration(1000)
		.attrTween("d", function(d) {
			this._current = this._current || d;
			var interpolate = d3.interpolate(this._current, d);
			this._current = interpolate(0);
			return function(t) {
				return donutData.arc(interpolate(t));
			};
		});

    slice.exit()
        .remove();

	function mousemoveHandler() {
            donutData.div.style("left", d3.event.pageX + 10 + "px");
            donutData.div.style("top", d3.event.pageY - 25 + "px");
            donutData.div.style("display", "inline-block");
            donutData.div.style("opacity", 0.8);
            donutData.div.html((this.__data__.data.label) + "<br>" + percentFormat(this.__data__.data.value) );

    }




	/* ------- TEXT LABELS -------*/

	var text = svg.select(".sliceLabels").selectAll("text")
		.data(donutData.pie(data), donutData.key);
		
	function midAngle(d){
		return d.startAngle + (d.endAngle - d.startAngle)/2;
	}

	text.enter()
		.append("text")
		.attr("dy", ".35em")
		.text(function(d) {
			return d.data.label;
		})
	  .merge(text)
	  .transition().duration(1000)
		.attrTween("transform", function(d) {
			this._current = this._current || d;
			var interpolate = d3.interpolate(this._current, d);
			this._current = interpolate(0);
			return function(t) {
				var d2 = interpolate(t);
				var pos = donutData.outerArc.centroid(d2);
				pos[0] = donutData.radius * (midAngle(d2) < Math.PI ? 1 : -1);
				return "translate("+ pos +")";
			};
		})
		.styleTween("text-anchor", function(d){
			this._current = this._current || d;
			var interpolate = d3.interpolate(this._current, d);
			this._current = interpolate(0);
			return function(t) {
				var d2 = interpolate(t);
				return midAngle(d2) < Math.PI ? "start":"end";
			};
		});

	text.exit()
		.remove();

	/* ------- SLICE TO TEXT POLYLINES -------*/

	var polyline = svg.select(".lines").selectAll("polyline")
		.data(donutData.pie(data), donutData.key);
	
	polyline.enter()
		.append("polyline")
		.merge(polyline)
    .transition().duration(1000)
		.attrTween("points", function(d){
			this._current = this._current || d;
			var interpolate = d3.interpolate(this._current, d);
			this._current = interpolate(0);
			return function(t) {
				var d2 = interpolate(t);
				var pos = donutData.outerArc.centroid(d2);
				pos[0] = donutData.radius * 0.95 * (midAngle(d2) < Math.PI ? 1 : -1);
				return [donutData.arc.centroid(d2), donutData.outerArc.centroid(d2), pos];
			};			
		});
	
	polyline.exit()
		.remove();
    
    
    // Toolip
    
    // function that creates and adds the tool tip to a selected element
    function tooltip(selection) {

        // add tooltip (svg circle element) when mouse enters label or slice
        selection.on('mouseenter', function (d) {

            donutData.svg.append('text')
                .attr('class', 'toolCircle')
                .attr('dy', -5) // hard-coded. can adjust this to adjust text vertical alignment in tooltip
                .html(toolTipHTML(d)) // add text to the circle.
                .style('font-size', '.7em')
                .style('text-anchor', 'middle'); // centres text in tooltip

            donutData.svg.append('circle')
                .attr('class', 'toolCircle')
                .attr('r', donutData.radius * 0.38) // radius of tooltip circle
                .style('fill', donutData.colorMap.get(d.data.label).color) // colour based on category mouse is over
                .style('fill-opacity', 0.65);

        });

        function toolTipHTML(d) {

            var tip = '',
                i   = 0;

            for (var key in d.data) {

                // if value is a number, format it as a percentage
                var value = (!isNaN(parseFloat(d.data[key]))) ? percentFormat(d.data[key]) : d.data[key];

                // leave off 'dy' attr for first tspan so the 'dy' attr on text element works. The 'dy' attr on
                // tspan effectively imitates a line break.
                if (i === 0) tip += '<tspan x="0">' + key + ': ' + value + '</tspan>';
                else tip += '<tspan x="0" dy="1.2em">' + key + ': ' + value + '</tspan>';
                i++;
            }

            return tip;
        }

        // color(d.data.label)
        // remove the tooltip when mouse leaves the slice/label
        selection.on('mouseout', function () {
            d3.selectAll('.toolCircle').remove();
            donutData.div.style("display", "none");
        });
    }

    d3.selectAll('path.slice').call(tooltip);



	// relax the label! Taken from https://codepen.io/anon/pen/jvpdeP. Similar to https://jsfiddle.net/thudfactor/HdwTH/
	// Needs more work,
    var alpha = 0.5,
        spacing = 15;

    function relax() {
        var again = false;
        text.each(function(d, i) {
            var a = this,
                da = d3.select(a),
                y1 = da.attr('y');
            text.each(function(d, j) {
                var b = this;
                if (a === b) {
                    return ;
                }

                db = d3.select(b);
                if (da.attr('text-anchor') !== db.attr('text-anchor')) {
                    return ;
                }

                var y2 = db.attr('y');
                deltaY = y1 - y2;

                if (Math.abs(deltaY) > spacing) {
                    return ;
                }

                again = true;
                sign = deltaY > 0? 1: -1;
                var adjust = sign * alpha;
                da.attr('y', +y1 + adjust);
                db.attr('y', +y2 - adjust);
            });
        });

        if (again) {
            var labelElements = text[0];
            polyline.attr('y2', function(d, i) {
                var labelForLine = d3.select(labelElements[i]);
                return labelForLine.attr('y');
            });
            setTimeout(relax, 20);
        }
    }

    //relax();


	return svg

};