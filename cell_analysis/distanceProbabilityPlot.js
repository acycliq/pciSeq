export class DistanceProbabilityPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;  // Store the shared tooltip
        this.margin = {top: 60, right: 80, bottom: 50, left: 100};
        this.initializePlot();
    }

    initializePlot() {
        // Get container width and height
        const container = d3.select(`#${this.containerId}`);
        const width = window.innerWidth - this.margin.left - this.margin.right;
        const height = window.innerHeight * 0.34 - this.margin.top - this.margin.bottom;

        // Create SVG
        this.svg = container.append('svg')
            .attr('width', width + this.margin.left + this.margin.right)
            .attr('height', height + this.margin.top + this.margin.bottom)
            .append('g')
            .attr('transform', `translate(${this.margin.left},${this.margin.top})`);

        // Set up scales
        this.x = d3.scaleLinear()
            .domain(d3.extent(this.data.x))
            .range([0, width]);

        this.y = d3.scaleLinear()
            .domain(d3.extent(this.data.y))
            .range([height, 0]);

        // Add axes
        this.svg.append('g')
            .attr('transform', `translate(0,${height})`)
            .call(d3.axisBottom(this.x));

        this.svg.append('g')
            .call(d3.axisLeft(this.y));

        // Add title
        this.svg.append("text")
            .attr("class", "plot-title")
            .attr("x", width / 2)
            .attr("y", -this.margin.top / 2)
            .style("text-anchor", "middle")
            .style("font-size", "16px")
            .text(`Distance vs Assignment Probability: Cell ${this.data.cell_num}`);

        // Add labels
        this.svg.append("text")
            .attr("x", width / 2)
            .attr("y", height + this.margin.bottom - 10)
            .style("text-anchor", "middle")
            .text(`Distance from cell ${this.data.cell_num} centroid`);

        this.svg.append("text")
            .attr("transform", "rotate(-90)")
            .attr("x", -height / 2)
            .attr("y", -this.margin.left + 40)
            .style("text-anchor", "middle")
            .text(`Assignment probability to cell ${this.data.cell_num}`);

        // // Create tooltip
        // this.tooltip = d3.select("body").append("div")
        //     .attr("class", "tooltip")
        //     .style("opacity", 0);

        this.updatePlot();
    }

    updatePlot() {
        // Create points data
        const points = this.data.x.map((d, i) => ({
            x: d,
            y: this.data.y[i],
            label: this.data.labels[i]
        }));

        // Add points
        const dots = this.svg.selectAll('circle')
            .data(points);

        // Enter new points
        dots.enter()
            .append('circle')
            .attr('r', 4)
            .style('fill', '#69b3a2')
            .merge(dots)
            .attr('cx', d => this.x(d.x))
            .attr('cy', d => this.y(d.y))
            .on('mouseover', (event, d) => {
                this.tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);
                this.tooltip.html(
                    `<strong>${d.label}</strong><br/>` +
                    `Distance: ${d.x.toFixed(2)}<br/>` +
                    `Probability: ${d.y.toFixed(2)}`
                )
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 28) + "px");
            })
            .on('mouseout', () => {
                this.tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            });

        // Remove old points
        dots.exit().remove();
    }

    resize() {
        // Implement resize logic here if needed
        // Similar to the existing ScatterPlot resize logic
    }
}
