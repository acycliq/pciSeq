import {calculateDimensions, PLOT_CONFIG} from "./plotConfig.js";
import { InterpretationGuide } from './interpretationGuide.js';

export class DistanceProbabilityPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;
        this.margin = {top: 60, right: 80, bottom: 50, left: 100};
        this.guide = null;

        // Create radius scale for contributions
        this.radiusScale = d3.scaleSqrt()
            .domain([0, d3.max(data.assigned_contributions) || 1]) // fallback if no contributions
            .range([0.1, 0.5]);  // min radius 3px, max 15px

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

        // Add x-axis label
        this.svg.append("text")
            .attr("x", width / 2)
            .attr("y", height + this.margin.bottom - 10)
            .style("text-anchor", "middle")
            .text(`Distance from cell ${this.data.cell_num} centroid`);

        // Add y-axis label
        this.svg.append("text")
            .attr("transform", "rotate(-90)")
            .attr("x", -height / 2)
            .attr("y", -this.margin.left + 40)
            .style("text-anchor", "middle")
            .text(`Assignment probability to cell ${this.data.cell_num}`);

        // Add interpretation guide
        this.guide = new InterpretationGuide(this.svg, this.width, this.height);
        this.guide.update('', '');

        const guideText = [
            "Notes:",
            `• Shows all spots in neighborhood of cell ${this.data.cell_num}`,
            `• Points close to x=0: Spots near cell ${this.data.cell_num} centroid`,
            `• Points with high y-values: Spots likely assigned to cell ${this.data.cell_num}`,
            `• Points with low y-values: Spots likely assigned to neighboring cells`,
            `• Point size indicates gene contribution`
        ];

        let guide = this.svg.select('.interpretation-guide');
        guide.selectAll("text").remove();
        guide.selectAll("text")
            .data(guideText)
            .enter()
            .append("text")
            .attr("x", 0)
            .attr("y", (d, i) => i * this.guide.lineHeight)
            .style("text-anchor", "start")
            .style("font-size", "12px")
            .text(d => d);

        // Position guide
        const guideBBox = guide.node().getBBox();
        guide.attr("transform",
            `translate(${this.width - guideBBox.width - 20}, ${20})`);

        this.updatePlot();
    }

    updatePlot() {
        // Create points data with all necessary attributes
        const points = this.data.x.map((x, i) => ({
            x: x,
            y: this.data.y[i],
            label: this.data.labels[i],
            contribution: this.data.assigned_contributions[i] || 0
        }));

        // Update points using join pattern
        const dots = this.svg.selectAll("circle")
            .data(points)
            .join("circle")
            .attr("cx", d => this.x(d.x))
            .attr("cy", d => this.y(d.y))
            .attr("r", d => this.radiusScale(d.contribution))
            .attr("fill", PLOT_CONFIG.point.color)
            .attr("stroke", "white")     // Add white outline
            .attr("stroke-width", "0.5") // Thin outline
            .on("mouseenter", (event, d) => {
                d3.select(event.target)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr("r", d => this.radiusScale(d.contribution) * 1.6)
                    .attr("fill", PLOT_CONFIG.point.hoverColor);
            })
            .on("mouseleave", (event, d) => {
                d3.select(event.target)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr("r", d => this.radiusScale(d.contribution))
                    .attr("fill", PLOT_CONFIG.point.color);
            })
            .on("mouseover", (event, d) => {
                this.tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);
                this.tooltip.html(
                    `<strong>${d.label}</strong><br/>` +
                    `Distance: ${d.x.toFixed(2)}<br/>` +
                    `Probability: ${d.y.toFixed(2)}<br/>` +
                    `Contribution: ${d.contribution.toFixed(2)}`
                )
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 28) + "px");
            })
            .on("mouseout", () => {
                this.tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            });
    }

    resize() {
        const { width, height } = calculateDimensions();

        // Update dimensions
        this.width = width;
        this.height = height;

        // Update SVG size
        this.svg.parent()
            .attr('width', width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom);

        // Update scales
        this.x.range([0, width]);
        this.y.range([height, 0]);

        // Update axes
        this.svg.select('.x-axis')
            .attr('transform', `translate(0,${height})`)
            .call(d3.axisBottom(this.x));

        this.svg.select('.y-axis')
            .call(d3.axisLeft(this.y));

        // Update points
        this.updatePlot();
    }
}
