import { PLOT_CONFIG, calculateDimensions } from './plotConfig.js';
import { InterpretationGuide } from './interpretationGuide.js';

export class ClassComparisonPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;
        this.defaultRadius = 5.5;
        this.visibleGenes = new Set(data.gene_names);
        this.currentClass = data.assigned_class;
        this.showRegression = false;  // Add this line

        this.currentAssignedClass = data.assigned_class;  // Fixed x-axis class
        this.currentUserClass = this.currentAssignedClass;      // Selected y-axis class

        this.setup();
        this.setupClassSelector();  // Single selector
        this.updatePlot(true);
    }

    createSvg() {
        return d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', this.width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', this.height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom)
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top})`);
    }

    createScales(width, height) {
        // Get index of current assigned class (for x-axis)
        const classIndex = this.data.class_names.indexOf(this.currentAssignedClass);

        // Get all x values (single cell reference)
        const xValues = this.data.scRNAseq_gene_counts.map(row => row[classIndex]);
        const xValues_adj = xValues.map((value, index) => value * this.data.gene_efficiency[index]);
        const xExtent = d3.extent(xValues_adj);

        // Get all y values (estimated counts for current selected class)
        const selectedClassIndex = this.data.class_names.indexOf(this.currentUserClass);
        const yValues = this.data.estimated_class_gene_counts_all.map(row => row[selectedClassIndex]);
        const yExtent = d3.extent(yValues);

        return {
            x: d3.scaleLinear().domain(xExtent).range([0, width]),
            y: d3.scaleLinear().domain(yExtent).range([height, 0])
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
        this.title = this.svg.append('text')
            .attr('class', 'title')
            .attr('text-anchor', 'middle')
            .attr('x', this.width / 2)
            .attr('y', -30)
            .style('font-size', '14px');

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

        this.updateLabels();
    }

    updateLabels() {
        this.title.text(`Class gene expressions`);

        // Update x-axis label with assigned class
        this.xLabel.text(`${this.data.assigned_class}: Single cell data`);

        // Update y-axis label with current selected class
        this.yLabel.text(`${this.currentClass}: Estimated gene counts`);

        this.subtitle.text('Comparing estimated vs scRNAseq gene counts');
    }

    setupDiagonalLine() {
        this.diagonalLine = this.svg.append('line')
            .attr('class', 'diagonal-line')
            .style('stroke', 'gray')
            .style('stroke-width', '1px')
            .style('stroke-dasharray', '4');

        this.updateDiagonalLine();
    }

    calculateRegression(plotData) {
        // Extract x and y values
        const x = plotData.map(d => d.x_adj);
        const y = plotData.map(d => d.y);

        // Calculate means
        const xMean = d3.mean(x);
        const yMean = d3.mean(y);

        // Calculate coefficients
        let ssxx = 0, ssyy = 0, ssxy = 0;
        for (let i = 0; i < x.length; i++) {
            ssxx += (x[i] - xMean) * (x[i] - xMean);
            ssyy += (y[i] - yMean) * (y[i] - yMean);
            ssxy += (x[i] - xMean) * (y[i] - yMean);
        }

        const slope = ssxy / ssxx;
        const intercept = yMean - slope * xMean;

        // Calculate R²
        const rSquared = (ssxy * ssxy) / (ssxx * ssyy);

        return { slope, intercept, rSquared };
    }

    setupRegressionControls() {
        const checkboxGroup = this.svg.append("g")
            .attr("class", "checkbox-group")
            .attr("transform", `translate(${this.width - 250}, -25)`);

        checkboxGroup.append("rect")
            .attr("width", 240)
            .attr("height", 30)
            .attr("fill", "white")
            .attr("opacity", 0.8)
            .attr("rx", 4);

        const foreignObject = checkboxGroup.append("foreignObject")
            .attr("width", 240)
            .attr("height", 30);

        const div = foreignObject.append("xhtml:div")
            .style("font-size", "12px")
            .style("color", "#666")
            .style("padding", "6px")
            .style("display", "flex")
            .style("align-items", "center")
            .style("gap", "10px");

        const label = div.append("label")
            .style("display", "flex")
            .style("align-items", "center")
            .style("gap", "5px");

        label.append("input")
            .attr("type", "checkbox")
            .on("change", (event) => {
                this.showRegression = event.target.checked;
                this.updatePlot();
            });

        label.append("span")
            .text("Draw regression line");

        // Add R² text next to checkbox label
        this.rSquaredText = div.append("span")
            .style("opacity", 0)
            .style("margin-left", "0px");
    }

    setupClassSelector() {
        // Create dropdown container
        const dropdownGroup = this.svg.append("g")
            .attr("class", "dropdown-group")
            .attr("transform", `translate(10, -45)`);

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

        // Add select element
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

        // Add options
        select.selectAll("option")
            .data(classOptions)
            .enter()
            .append("option")
            .attr("value", d => d.name)
            .property("selected", d => d.name === this.currentUserClass)
            .text(d => d.text);

        // Add change handler
        select.on("change", (event) => {
            this.currentUserClass = event.target.value;  // Update the correct variable
            const classIndex = this.data.class_names.indexOf(this.currentUserClass);
            const yValues = this.data.estimated_class_gene_counts_all.map(d => d[classIndex]);

            // Only update y scale with new class data
            const yExtent = d3.extent(yValues);
            this.scales.y.domain(yExtent);

            // Update only y-axis with animation
            this.yAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisLeft(this.scales.y));

            this.updatePlot(false);
        });
    }

    setup() {
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        this.svg = this.createSvg();
        this.scales = this.createScales(width, height);

        this.guide = new InterpretationGuide(this.svg, width, height, true);

        this.setupAxes();
        this.setupLabels();
        this.setupDiagonalLine();
        this.setupRegressionLine();
        this.setupRegressionControls();  // Add this line
    }

    setupRegressionLine() {
        // Add regression line with soft red color and dotted style
        this.regressionLine = this.svg.append('line')
            .attr('class', 'regression-line')
            .style('stroke', '#FF9999')  // Changed to soft red
            .style('stroke-width', '2px')
            .style('stroke-dasharray', '3,3')  // Added dotted line style
            .style('opacity', 0);  // Start hidden

        // Add R² text
        this.rSquaredText = this.svg.append('text')
            .attr('class', 'r-squared')
            .attr('text-anchor', 'end')
            .attr('x', this.width - 10)
            .attr('y', 30)
            .style('font-size', '12px')
            .style('opacity', 0);  // Start hidden
    }


    updateRegressionLine(plotData) {
        if (!this.showRegression) {
            this.regressionLine.style('opacity', 0);
            this.rSquaredText.style('opacity', 0);
            return;
        }

        const { slope, intercept, rSquared } = this.calculateRegression(plotData);

        const xDomain = this.scales.x.domain();

        this.regressionLine
            .style('opacity', 1)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('x1', this.scales.x(xDomain[0]))
            .attr('y1', this.scales.y(slope * xDomain[0] + intercept))
            .attr('x2', this.scales.x(xDomain[1]))
            .attr('y2', this.scales.y(slope * xDomain[1] + intercept));

        // Update R² text in the checkbox group
        this.rSquaredText
            .style('opacity', 1)
            .text(`R² = ${rSquared.toFixed(3)}`);
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

        dots.merge(dotsEnter)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('cx', d => this.scales.x(d.x_adj))
            .attr('cy', d => this.scales.y(d.y))
            .attr('r', this.defaultRadius)
            .attr('fill', PLOT_CONFIG.point.color);

        this.svg.selectAll('circle')
            .on('mouseenter', (event, d) => {
                const circle = d3.select(event.target);
                circle.transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr('r', this.defaultRadius * 1.6)
                    .attr('fill', PLOT_CONFIG.point.hoverColor);
            })
            .on('mouseleave', (event, d) => {
                const circle = d3.select(event.target);
                circle.transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr('r', this.defaultRadius)
                    .attr('fill', PLOT_CONFIG.point.color);
            })
            .on('mouseover', (event, d) => {
                const i = this.data.gene_names.indexOf(d.name);
                const efficiency = this.data.gene_efficiency[i];
                const vw = Math.max(document.documentElement.clientWidth || 0, window.innerWidth || 0);

                this.tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);

                this.tooltip.html(
                    `<strong>${d.name}</strong><br>` +
                    `scRNAseq (adjusted): ${d.x_adj.toFixed(3)}<br>` +
                    `scRNAseq (raw): ${d.x_adj.toFixed(3)}<br>` +
                    `Observed: ${d.y.toFixed(3)}<br>` +
                    `Efficiency: ${efficiency.toFixed(3)}`
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
                    .style("top", `${event.pageY - 28}px`);
            })
            .on("mouseout", () => {
                this.tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            });
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

    updateClass(newClass) {
        this.currentClass = newClass;
        this.updatePlot(false);
    }

    updatePlot(initial = false) {
        // Get indices for both classes
        const assignedClassIndex = this.data.class_names.indexOf(this.currentAssignedClass);
        const selectedClassIndex = this.data.class_names.indexOf(this.currentUserClass);

        const plotData = this.data.gene_names
            .map((name, i) => ({
                name,
                x_adj: this.data.scRNAseq_gene_counts[i][assignedClassIndex] * this.data.gene_efficiency[i],  // x uses assigned class
                y: this.data.estimated_class_gene_counts_all[i][selectedClassIndex],  // y uses selected class
                g: this.data.gene_efficiency[i]
            }))
            .filter(d => this.visibleGenes.has(d.name));

        this.updatePoints(plotData);
        this.updateLabels();
        this.updateDiagonalLine();
        this.guide.update(this.scales);
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
