import {calculateDimensions, PLOT_CONFIG} from "./plotConfig.js";
import { InterpretationGuide } from './interpretationGuide.js';

export class DistanceProbabilityPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;
        this.margin = {top: 60, right: 80, bottom: 50, left: 100};
        this.guide = null;
        this.sizeByCount = false;
        this.defaultRadius = 5.5;
        this.visibleGenes = new Set(data.labels);

        this.radiusScale = d3.scaleSqrt()
            .domain([0, d3.max(Object.values(data.gene_counts)) || 1])
            .range([2, 8]);

        this.initializePlot();
    }

    initializePlot() {
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        this.svg = d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom)
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top})`);

        this.x = d3.scaleLinear()
            .domain(d3.extent(this.data.x))
            .range([0, width]);

        this.y = d3.scaleLinear()
            .domain(d3.extent(this.data.y))
            .range([height, 0]);

        this.setupAxes();
        this.setupLabels();
        this.setupGuide();
        this.setupSizeControl();
        this.updatePlot();
    }

    setupAxes() {
        this.svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${this.height})`)
            .call(d3.axisBottom(this.x)
                .ticks(5)
                .tickSize(-5)
                .tickPadding(5));

        this.svg.append('g')
            .attr('class', 'y-axis')
            .call(d3.axisLeft(this.y)
                .ticks(5)
                .tickSize(-5)
                .tickPadding(5));
    }

    setupLabels() {
        this.svg.append("text")
            .attr("class", "plot-title")
            .attr("x", this.width / 2)
            .attr("y", -this.margin.top / 2)
            .style("text-anchor", "middle")
            .style("font-size", "16px")
            .text(this.data.title);

        this.svg.append("text")
            .attr("class", "axis-label x-axis-label")
            .attr("x", this.width / 2)
            .attr("y", this.height + this.margin.bottom - 10)
            .style("text-anchor", "middle")
            .text(this.data.xlabel);

        this.svg.append("text")
            .attr("class", "axis-label y-axis-label")
            .attr("transform", "rotate(-90)")
            .attr("x", -this.height / 2)
            .attr("y", -this.margin.left + 40)
            .style("text-anchor", "middle")
            .text(this.data.ylabel);
    }

    setupGuide() {
        this.guide = new InterpretationGuide(this.svg, this.width, this.height);
        this.guide.update('', '');

        const guideText = [
            `• Shows all spots in neighborhood of cell ${this.data.cell_num}`,
            `• Points close to x=0: Spots near cell ${this.data.cell_num} centroid`,
            `• Points with high y-values: Spots likely assigned to cell ${this.data.cell_num}`,
            `• Points with low y-values: Spots likely assigned to neighboring cells`
        ];

        let guide = this.svg.select('.interpretation-guide');
        guide.selectAll("text").remove();
        guide.selectAll("text")
            .data(guideText)
            .enter()
            .append("text")
            .attr("x", 0)
            .attr("y", (d, i) => i * 15)
            .style("text-anchor", "start")
            .style("font-size", "12px")
            .text(d => d);

        const guideBBox = guide.node().getBBox();
        guide.attr("transform", `translate(${this.width - guideBBox.width - 20}, 40)`);
    }

    setupSizeControl() {
        const checkboxGroup = this.svg.append("g")
            .attr("class", "checkbox-group")
            .attr("transform", `translate(${this.width - 200}, -5)`);

        checkboxGroup.append("rect")
            .attr("width", 180)
            .attr("height", 20)
            .attr("fill", "white")
            .attr("opacity", 0.8)
            .attr("rx", 4);

        checkboxGroup.append("foreignObject")
            .attr("width", 180)
            .attr("height", 20)
            .append("xhtml:div")
            .style("font-size", "12px")
            .html(`
                <input type="checkbox" id="size-toggle-${this.containerId}">
                <label for="size-toggle-${this.containerId}">Size by gene count</label>
            `);

        d3.select(`#size-toggle-${this.containerId}`)
            .on("change", (event) => {
                this.sizeByCount = event.target.checked;
                this.updatePlot();
            });
    }

    updateVisibleGenes(visibleGenes) {
        this.visibleGenes = new Set(visibleGenes);
        this.updatePlot();
    }

    updatePlot() {
        const points = this.data.x.map((x, i) => ({
            x: x,
            y: this.data.y[i],
            label: this.data.labels[i],
            geneCount: this.data.gene_counts[this.data.labels[i]] || 0
        })).filter(d => this.visibleGenes.has(d.label));

        const dots = this.svg.selectAll("circle")
            .data(points, d => d.label);

        dots.exit()
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr("r", 0)
            .remove();

        const dotsEnter = dots.enter()
            .append("circle")
            .attr("fill", PLOT_CONFIG.point.color)
            .attr("stroke", "white")
            .attr("stroke-width", "0.5")
            .attr("cx", d => this.x(d.x))
            .attr("cy", d => this.y(d.y))
            .attr("r", 0);

        dots.merge(dotsEnter)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr("cx", d => this.x(d.x))
            .attr("cy", d => this.y(d.y))
            .attr("r", d => this.sizeByCount ?
                this.radiusScale(d.geneCount) :
                this.defaultRadius);

        this.svg.selectAll("circle")
            .on("mouseenter", (event, d) => {
                d3.select(event.target)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr("r", d => this.sizeByCount ?
                        this.radiusScale(d.geneCount) * 1.6 :
                        this.defaultRadius * 1.6)
                    .attr("fill", PLOT_CONFIG.point.hoverColor);
            })
            .on("mouseleave", (event, d) => {
                d3.select(event.target)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr("r", d => this.sizeByCount ?
                        this.radiusScale(d.geneCount) :
                        this.defaultRadius)
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
                    `Gene Count: ${d.geneCount.toFixed(2)}`
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
        this.width = width;
        this.height = height;

        const svg = d3.select(`#${this.containerId}`).select('svg')
            .attr('width', width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom);

        this.x.range([0, width]);
        this.y.range([height, 0]);

        this.svg.select('.x-axis')
            .attr('transform', `translate(0,${height})`)
            .call(d3.axisBottom(this.x));

        this.svg.select('.y-axis')
            .call(d3.axisLeft(this.y));

        this.updatePlot();
    }
}
