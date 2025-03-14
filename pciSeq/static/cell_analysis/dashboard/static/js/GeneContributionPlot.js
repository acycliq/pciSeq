import { PLOT_CONFIG, calculateDimensions } from './plotConfig.js';
import { InterpretationGuide } from './interpretationGuide.js';

export class GeneContributionPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.currentUserClass = data.user_class;
        this.currentAssignedClass = data.assigned_class;
        this.tooltip = tooltip;
        this.sizeByCount = false;
        this.defaultRadius = 5.5;
        this.visibleGenes = new Set(data.gene_names);

        // NEW: Add ratio coloring properties
        this.colorByRatio = false;
        this.ratioColorScale = d3.scaleLinear()
            .domain([0, 1, 2])
            .range(['#2166AC', '#4D4D4D', '#B2182B']);  // Blue - Gray - Red
            // .range(['blue', 'black', 'red']);

        this.radiusScale = d3.scaleSqrt()
            .domain([0, d3.max(data.gene_counts) || 1])
            .range([2, 8]);

        this.setup();
        this.setupClassSelector();
        this.updatePlot(true);
    }

    // NEW: Add ratio calculation method
    calculateRatio(geneCount, ClassAvgCounts_pciSeq = 1) {
        if (ClassAvgCounts_pciSeq === 0) return 0;
        // const classGeneCountAdj = classGeneCount * geneInefficiency;
        return geneCount / ClassAvgCounts_pciSeq;
    }

    createSvg() {
        return d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', this.width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', this.height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom + 50)  // Added 40px for dropdown
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top + 50})`);  // Added 40 to top margin
    }

    createScales(width, height) {
        const xContributions = this.data.contr[this.currentAssignedClass];
        const xExtent = d3.extent(xContributions);

        const yContributions = this.data.contr[this.currentUserClass];
        const yExtent = d3.extent(yContributions);

        return {
            x: d3.scaleLinear().domain(xExtent).range([0, width]),   // x scale uses only assigned class
            y: d3.scaleLinear().domain(yExtent).range([height, 0])   // y scale uses only selected class
        };
    }

    setupAxes() {
        this.xAxis = this.svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${this.height})`)
            .call(d3.axisBottom(this.scales.x));

        this.yAxis = this.svg.append('g')
            .attr('class', 'y-axis')
            .call(d3.axisLeft(this.scales.y));
    }

    setupLabels() {
        this.xLabel = this.svg.append('text')
            .attr('class', 'x-label')
            .attr('text-anchor', 'middle')
            .attr('x', this.width / 2)
            .attr('y', this.height + PLOT_CONFIG.margin.bottom - 5);

        this.yLabel = this.svg.append('text')
            .attr('class', 'y-label')
            .attr('text-anchor', 'middle')
            .attr('transform', 'rotate(-90)')
            .attr('x', -this.height / 2)
            .attr('y', -PLOT_CONFIG.margin.left + 15);

        this.subtitle = this.svg.append('text')
            .attr('class', 'subtitle')
            .attr('text-anchor', 'middle')
            .attr('x', this.width / 2)
            .attr('y', -10)
            .style('font-size', '12px');

        this.title = this.svg.append('text')
            .attr('class', 'title')
            .attr('text-anchor', 'middle')
            .attr('x', this.width / 2)
            .attr('y', -30)
            .style('font-size', '14px')
            .text(`Cell ${this.data.cell_num}: Gene likelihood contributions`);

        this.updateLabels();
    }

    setupDiagonalLine() {
        this.diagonalLine = this.svg.append('line')
            .attr('class', 'diagonal-line')
            .style('stroke', 'gray')
            .style('stroke-width', '1px')
            .style('stroke-dasharray', '4');

        this.updateDiagonalLine();
    }

    setupSizeControl() {
        const checkboxGroup = this.svg.append("g")
            .attr("class", "checkbox-group")
            .attr("transform", `translate(${this.width - 200}, -110)`); //from -25

        checkboxGroup.append("rect")
            .attr("width", 180)
            .attr("height", 40)  // Increased height for two checkboxes
            .attr("fill", "white")
            .attr("opacity", 0.8)
            .attr("rx", 4);

        // Original checkbox (hidden)
        checkboxGroup.append("foreignObject")
            .attr("width", 180)
            .attr("height", 20)
            .append("xhtml:div")
            .style("font-size", "12px")
            .style("display", "none")  // Hidden
            .html(`
                <input type="checkbox" id="size-toggle-${this.containerId}">
                <label for="size-toggle-${this.containerId}">Size by gene count</label>
            `);

        // NEW: Add ratio coloring checkbox
        checkboxGroup.append("foreignObject")
            .attr("width", 180)
            .attr("height", 20)
            .append("xhtml:div")
            .style("font-size", "12px")
            .html(`
                <input type="checkbox" id="ratio-color-toggle-${this.containerId}">
                <label for="ratio-color-toggle-${this.containerId}">Color by expression ratio</label>
            `);

        // Original handler (kept but hidden)
        d3.select(`#size-toggle-${this.containerId}`)
            .on("change", (event) => {
                this.sizeByCount = event.target.checked;
                this.updatePlot(false);
            });

        // NEW: Add ratio coloring handler
        d3.select(`#ratio-color-toggle-${this.containerId}`)
            .on("change", (event) => {
                this.colorByRatio = event.target.checked;
                this.colorLegend.transition()
                    .duration(PLOT_CONFIG.animation.duration)
                    .style('opacity', this.colorByRatio ? 1 : 0);
                this.updatePlot(false);
            });
    }

    // NEW: Add color legend setup
    setupColorLegend() {
        const legendWidth = 300;
        const legendHeight = 40;
        const margin = 10;

        const legend = this.svg.append('g')
            .attr('class', 'color-legend')
            .attr('transform', `translate(${this.width - legendWidth - margin}, ${-90})`)
            .style('opacity', 0);

        const gradient = legend.append('defs')
            .append('linearGradient')
            .attr('id', `ratio-gradient-${this.containerId}`)
            .attr('x1', '0%')
            .attr('x2', '100%');

        gradient.append('stop')
            .attr('offset', '0%')
            .attr('stop-color', '#2166AC');
        gradient.append('stop')
            .attr('offset', '50%')
            .attr('stop-color', '#4D4D4D');
        gradient.append('stop')
            .attr('offset', '100%')
            .attr('stop-color', '#B2182B');

        legend.append('rect')
            .attr('width', legendWidth)
            .attr('height', legendHeight)
            .attr('fill', 'white')
            .attr('opacity', 0.8)
            .attr('rx', 4);

        legend.append('rect')
            .attr('x', 10)
            .attr('y', 10)
            .attr('width', legendWidth - 20)
            .attr('height', 10)
            .attr('fill', `url(#ratio-gradient-${this.containerId})`);

        legend.append('text')
            .attr('x', 10)
            .attr('y', 35)
            .attr('text-anchor', 'start')
            .style('font-size', '10px')
            .text('Under-expressed');

        legend.append('text')
            .attr('x', legendWidth/2)
            .attr('y', 35)
            .attr('text-anchor', 'middle')
            .style('font-size', '10px')
            .text('Expected');

        legend.append('text')
            .attr('x', legendWidth - 10)
            .attr('y', 35)
            .attr('text-anchor', 'end')
            .style('font-size', '10px')
            .text('Over-expressed');

        this.colorLegend = legend;
    }

    setup() {
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        this.svg = this.createSvg();
        this.scales = this.createScales(
            width,
            height,
            this.data.contr[this.currentAssignedClass],
            this.data.contr[this.currentUserClass]
        );

        this.guide = new InterpretationGuide(this.svg, width, height);

        this.setupAxes();
        this.setupLabels();
        this.setupDiagonalLine();
        this.setupSizeControl();
        this.setupColorLegend();  // NEW: Add color legend setup
    }

    setupClassSelector() {
        // Move dropdown further up to avoid title overlap
        const dropdownGroup = this.svg.append("g")
            .attr("class", "dropdown-group")
            .attr("transform", `translate(0, -85)`);  // Changed from -45 to -85

        // Add background rectangle
        dropdownGroup.append("rect")
            .attr("width", 200)
            .attr("height", 30)
            .attr("fill", "white")
            .attr("opacity", 0.8)
            .attr("rx", 4);

        // Create array of class options
        const classOptions = this.data.class_names.map(className => ({
            name: className,
            probability: this.data.class_probs[className],
            text: `${className} (${(this.data.class_probs[className] * 100).toFixed(2)}%)`
        })).sort((a, b) => b.probability - a.probability);

        // Add select element using foreignObject
        const select = dropdownGroup.append("foreignObject")
            .attr("width", 190)
            .attr("height", 30)
            .attr("x", 5)
            .attr("y", 2)
            .append("xhtml:select")
            .style("width", "100%")
            .style("height", "26px")
            .style("font-size", "12px")
            .style("border-radius", "4px")
            .style("border", "1px solid #ccc");

        // Add options to select
        select.selectAll("option")
            .data(classOptions)
            .enter()
            .append("option")
            .attr("value", d => d.name)
            .property("selected", d => d.name === this.currentUserClass)
            .text(d => d.text);

        // Add change handler
        select.on("change", (event) => {
            this.currentUserClass = event.target.value;

            // Recalculate scales with new class data
            this.scales = this.createScales(this.width, this.height);

            // Update y-axis with animation
            this.yAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisLeft(this.scales.y));

            // Update the rest of the plot
            this.updatePlot(false);
        });
    }

    updatePoints(plotData) {
        const dots = this.svg.selectAll('circle')
            .data(plotData, d => d.name);

        dots.exit()
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('r', 0)
            .remove();

        const dotsEnter = dots.enter()
            .append('circle')
            .attr('fill', PLOT_CONFIG.point.color)
            .attr('stroke', 'white')
            .attr('stroke-width', '0.5')
            .attr('r', 0);

        // MODIFIED: Update points with ratio coloring
        dots.merge(dotsEnter)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('cx', d => this.scales.x(d.x))
            .attr('cy', d => this.scales.y(d.y))
            .attr('r', d => this.sizeByCount ?
                this.radiusScale(d.geneCount) :
                this.defaultRadius)
            .attr('fill', d => {
                if (this.colorByRatio) {
                    const index = this.data.gene_names.indexOf(d.name);
                    const ratio = this.calculateRatio(
                        d.geneCount,
                        d.ClassAvgCounts_pciSeq
                    );
                    return this.ratioColorScale(ratio);
                }
                return PLOT_CONFIG.point.color;
            });

        // MODIFIED: Update hover behavior
        this.svg.selectAll('circle')
            .on('mouseenter', (event, d) => {
                const circle = d3.select(event.target);
                const currentFill = circle.attr('fill');

                circle.transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr('r', d => this.sizeByCount ?
                        this.radiusScale(d.geneCount) * 1.6 :
                        this.defaultRadius * 1.6)
                    .attr('fill', this.colorByRatio ?
                        currentFill :
                        PLOT_CONFIG.point.hoverColor);
            })
            .on('mouseleave', (event, d) => {
                const circle = d3.select(event.target);

                circle.transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr('r', d => this.sizeByCount ?
                        this.radiusScale(d.geneCount) :
                        this.defaultRadius)
                    .attr('fill', this.colorByRatio ?
                        this.ratioColorScale(this.calculateRatio(
                            d.geneCount,
                            d.ClassAvgCounts_pciSeq
                        )) :
                        PLOT_CONFIG.point.color);
            })
            .on('mouseover', (event, d) => {
                const vw = Math.max(document.documentElement.clientWidth || 0, window.innerWidth || 0);
                const ratio = this.calculateRatio(
                    d.geneCount,
                    d.ClassAvgCounts_pciSeq
                );

                this.tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);

                this.tooltip.html(
                    `<strong>${d.name}</strong><br>` +
                    `X: ${d.x.toFixed(3)}<br>` +
                    `Y: ${d.y.toFixed(3)}<br>` +
                    `Cell Gene Count: ${d.geneCount.toFixed(2)}<br>` +
                    `Class Avg Counts (pciSeq): ${d.ClassAvgCounts_pciSeq.toFixed(2)}<br>` +
                    `Class Avg Counts (scRNAseq): ${d.ClassAvgCounts_scRNAseq.toFixed(2)}<br>`
                    // `Expression Ratio: ${ratio.toFixed(2)}`
                );

                const tooltipWidth = this.tooltip.node().getBoundingClientRect().width;
                let left = event.pageX;

                if (event.pageX + tooltipWidth/2 > vw) {
                    left = vw - tooltipWidth - 10;
                }
                else if (event.pageX - tooltipWidth/2 < 0) {
                    left = tooltipWidth/2 + 10;
                }

                this.tooltip
                    .style("left", `${left}px`)
                    .style("top", `${event.pageY - 24}px`);
            })
            .on("mouseout", () => {
                this.tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            });
    }

    getAxisLabel(className, probability) {
        return `Likelihood contr to class: ${className} (${(probability * 100).toFixed(2)}%)`;
    }

    updateLabels() {
        this.yLabel.text(this.getAxisLabel(
            this.currentUserClass,
            this.data.class_probs[this.currentUserClass]
        ));

        // Keep x-label static since it's the assigned class
        this.xLabel.text(this.getAxisLabel(
            this.currentAssignedClass,
            this.data.class_probs[this.currentAssignedClass]
        ));

        this.subtitle.text(
            `Assigned class: ${this.currentAssignedClass} vs Selected class: ${this.currentUserClass}`
        );
    }

    updateDiagonalLine() {
        const minDomain = Math.min(
            this.scales.x.domain()[0],
            this.scales.y.domain()[0]
        );
        const maxDomain = Math.max(
            this.scales.x.domain()[1],
            this.scales.y.domain()[1]
        );

        this.diagonalLine.transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('x1', this.scales.x(minDomain))
            .attr('y1', this.scales.y(minDomain))
            .attr('x2', this.scales.x(maxDomain))
            .attr('y2', this.scales.y(maxDomain));
    }

    updateVisibleGenes(visibleGenes) {
        this.visibleGenes = visibleGenes;
        this.updatePlot(false);
    }

    updatePlot(initial = false) {
        const classIndex = this.data.class_names.indexOf(this.currentAssignedClass)
        const plotData = this.data.gene_names
            .map((name, i) => ({
                name,
                x: this.data.contr[this.currentAssignedClass][i],
                y: this.data.contr[this.currentUserClass][i],
                geneCount: this.data.gene_counts[i],
                ClassAvgCounts_pciSeq: this.data.mean_gene_reads_per_class[i][classIndex],
                ClassAvgCounts_scRNAseq: this.data.scRNAseq_gene_counts[i][classIndex]
            }))
            .filter(d => this.visibleGenes.has(d.name));

        // Remove scale recalculation for non-initial updates
        if (initial) {
            this.scales = this.createScales(this.width, this.height);

            this.xAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisBottom(this.scales.x));

            this.yAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisLeft(this.scales.y));
        }

        this.updatePoints(plotData);
        this.updateLabels();
        this.updateDiagonalLine();
        // this.guide.update(this.scales); // Hide the guide
    }

    resize() {
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        this.svg.attr('width', width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom);

        this.scales.x.range([0, width]);
        this.scales.y.range([height, 0]);

        this.updatePlot(false);
    }
}
