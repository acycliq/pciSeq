import {calculateDimensions, PLOT_CONFIG} from "./plotConfig.js";
import { InterpretationGuide } from './interpretationGuide.js';

export class DistanceProbabilityPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;  // Store the shared tooltip
        this.margin = {top: 60, right: 80, bottom: 50, left: 100};
        this.guide = null;
        this.initializePlot();
    }

    initializePlot() {
        // Get container width and height
        const container = d3.select(`#${this.containerId}`);

        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        // Create SVG
        this.svg = container.append('svg')
            .attr('width', width + this.margin.left + this.margin.right)
            .attr('width', this.width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', this.height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom)
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top})`);

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

        // Add interpretation guide with correct positioning
        this.guide = new InterpretationGuide(this.svg, this.width, this.height);
        // Pass empty strings as parameters since we're not using the class/assigned class paradigm
        this.guide.update('', '');

        // Now update with our custom text using the existing guide's updateGuideText method
        const guideText = [
            "Notes:",
            `• Shows all spots in neighborhood of cell ${this.data.cell_num}`,
            `• Points close to x=0: Spots near cell ${this.data.cell_num} centroid`,
            `• Points with high y-values: Spots likely assigned to cell ${this.data.cell_num}`,
            `• Points with low y-values: Spots likely assigned to neighboring cells`,
            // `• Pattern shows how cell territory influences spot assignments`
        ];

        let guide = this.svg.select('.interpretation-guide');
        guide.selectAll("text").remove();  // Remove existing text
        guide.selectAll("text")  // Add new text
            .data(guideText)
            .enter()
            .append("text")
            .attr("x", 0)
            .attr("y", (d, i) => i * this.guide.lineHeight)
            .style("text-anchor", "start")
            .style("font-size", "12px")
            .text(d => d);

        // Reposition guide to top-right corner
        const guideBBox = guide.node().getBBox();
        guide.attr("transform",
            `translate(${this.width - guideBBox.width - 20}, ${20})`);  // 20px margin from top and right


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
            .attr('r', PLOT_CONFIG.point.radius)
            .attr('fill', PLOT_CONFIG.point.color)
            .merge(dots)
            .attr('cx', d => this.x(d.x))
            .attr('cy', d => this.y(d.y))
            .on('mouseenter', function(event, d) {  // Note: need event parameter in D3v6+
                console.log('Mouse enter')
                d3.select(this)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr('r', 1.6 * PLOT_CONFIG.point.radius)
                    .attr('fill', '#4a90e2');
            })
            .on('mouseleave', function(event, d) {  // Note: need event parameter in D3v6+
                d3.select(this)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr('r', PLOT_CONFIG.point.radius)
                    .attr('fill', PLOT_CONFIG.point.color);
            })
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

    updateVisibleGenes(visibleGenes) {
    // Update points visibility
    this.svg.selectAll('circle')
        .transition()
        .duration(200)
        .style('opacity', d => visibleGenes.includes(d.name) ? 1 : 0.1)
        .style('pointer-events', d => visibleGenes.includes(d.name) ? 'all' : 'none');
    }

    resize() {
        // Implement resize logic here if needed
        // Similar to the existing ScatterPlot resize logic
    }
}
