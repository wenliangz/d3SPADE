/* fix overflow ( text outside the boundaries of the plot)
 * add marker selection and recoloring
 */

var width = 600,
    height = 600,
    margin = 10,
    link,
    node,
    colorBy='Reference',
    sizeBy;

var force = d3.layout.force()
    .size([width-margin, height-margin]);

var svg = d3.select("#chart").append("svg")
    .attr("width", width)
    .attr("height", height);

function draw(json) {
  // size by
  var sizeOptions = d3.select('#sizeBy');

  sizeOptions
    .selectAll('option')
      .data(json.stimulations)
      .enter()
      .append('option')
        .text(function(d) {return d;});

  sizeOptions
    .on('change',updateSize);

	// color by
	var colorOptions = d3.select('#colorBy');

	colorOptions
		.selectAll('option')
			.data(json.reagents)
			.enter()
			.append('option')
				.text(function(d) {return d;});

	colorOptions
		.on('change',updateColor);

  // define link weight scale
	var weight_extent = d3.extent(json.links,function(d) { return d.weight; });
	weight = d3.scale.linear()
	    .domain(weight_extent)
	    .range([1, 5]);

	// define x scale
	var x_extent = d3.extent(json.nodes,function(d) { return d.x_fixed; });

    var x_scale = d3.scale.linear()
        .range([0+margin,width-margin])
        .domain(x_extent);

    // define y scale
	var y_extent = d3.extent(json.nodes,function(d) { return d.y_fixed; });

    var y_scale = d3.scale.linear()
        .range([height-margin,0+margin])
        .domain(y_extent);

    json.nodes.forEach(function(d) {
		d.x = x_scale(d.x_fixed)
		d.y = y_scale(d.y_fixed)
		//d.fixed=1
	});

	force
	    .charge(-30)
	    .linkDistance(5)
		.nodes(json.nodes)
		.links(json.links)
		.start()
		.stop();

	link = svg.selectAll("line.link")
		.data(json.links)
		.enter()
		.append("line")
			.attr("class", "link")
			.attr("x1", function(d) { return d.source.x; })
			.attr("y1", function(d) { return d.source.y; })
			.attr("x2", function(d) { return d.target.x; })
			.attr("y2", function(d) { return d.target.y; })
			.style("stroke-width", function(d) { return weight(d.weight); });

	node = svg.selectAll("circle.node")
		.data(json.nodes)
		.enter()
		.append("circle")
			.attr("class", "node")
			.attr('cx', function(d) {return d.x })
			.attr('cy', function(d) {return d.y })
		.call(force.drag)
		.on('mouseover.tooltip',function(d) {
                d3.select('text#n'+d.name).remove();
                svg
                    .append('text')
                    .text(colorBy+': '+Math.round(100*d[sizeBy][colorBy])/100)
                    .attr('x',d.x)
                    .attr('y',d.y-10)
                    .attr('id','n'+d.name);
        })
        .on('mouseout.tooltip',function(d) {
            d3.select('text#n'+d.name)
                .remove();
        });

	force.on("tick", function() {
		link.attr("x1", function(d) { return d.source.x; })
			.attr("y1", function(d) { return d.source.y; })
			.attr("x2", function(d) { return d.target.x; })
			.attr("y2", function(d) { return d.target.y; });

		node.attr("cx", function(d) { return d.x }) //if(d.x > width-margin || d.x < margin) { return margin } else { return d.x }
			.attr("cy", function(d) { return d.y }); //if(d.y > height-margin || d.y < margin) { return margin} else { return d.y }
	});

  updateSize();

	updateColor();

  function updateSize() {

		sizeBy  = sizeOptions.node().options[sizeOptions.node().options.selectedIndex].value;

    // define color scale
		var size_values = json.nodes.map(function(d) { return Math.pow(d[sizeBy]['percenttotal'],0.57);}),
        size_extent = d3.extent(size_values);

		var size = d3.scale.linear()
		    	.domain(size_extent)
		    	.range([0,12]);

		node.attr("r", function(d) { return size(Math.pow(d[sizeBy]['percenttotal'],0.57)); });
    //console.log(sizeBy);
    updateColor();

		return sizeBy;

	}

	function updateColor() {

		colorBy  = colorOptions.node().options[colorOptions.node().options.selectedIndex].value;

    // define color scale
		var col_values = json.nodes.map(function(d) { return d[sizeBy][colorBy];}),
			col_extent = [];

    col_values = col_values.map(function(d) {if(!isNaN(d)) { return(d)}} )

		col_values.sort(function(a,b) { return a-b;});

		col_extent.push(d3.min(col_values));
		col_extent.push(d3.quantile(col_values,0.25));
		col_extent.push(d3.quantile(col_values,0.75));
		col_extent.push(d3.max(col_values));

		var color = d3.scale.linear()
		    	.domain(col_extent)
		    	.range(["blue", "yellow", "green", "red"]);

		node.style("fill", function(d) { if(d[sizeBy][colorBy]=='NA') { return('white') } else { return color(d[sizeBy][colorBy]); } });
    //console.log(colorBy);
    //console.log(col_extent);

		return colorBy;

	}

}
