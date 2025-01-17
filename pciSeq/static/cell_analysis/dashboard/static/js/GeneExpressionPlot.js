import { PLOT_CONFIG, calculateDimensions } from './plotConfig.js';  // or wherever it's defined

export class GeneExpressionPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;
        this.defaultRadius = 5.5;
        this.visibleGenes = new Set(data.gene_names);
        this.showRegression = false;
        this.currentUserClass = this.data.user_class;

        this.setup();
        this.updatePlot(true);
    }

    createSvg() {
        return d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', this.width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', this.height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom + 40)
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top + 40})`);
    }

    createScales(width, height) {
        // Get index of current class
        const classIndex = this.data.class_names.indexOf(this.currentUserClass);

        // Get the data for both axes
        const xValues = this.data.gene_counts;

        // Get all y values (single cell reference)
        const yValues = this.data.scRNAseq_gene_counts.map(row => row[classIndex]);
        const yValues_adj = yValues.map((value, index) => value * this.data.gene_efficiency[index]);

        // Calculate extents separately for x and y
        const xExtent = d3.extent(xValues);
        const yExtent = d3.extent(yValues_adj);

        // Add some padding to the extents
        const xPadding = (xExtent[1] - xExtent[0]) * 0.05;
        const yPadding = (yExtent[1] - yExtent[0]) * 0.05;

        return {
            x: d3.scaleLinear()
                .domain([xExtent[0], xExtent[1]])
                .range([0, width]),
            y: d3.scaleLinear()
                .domain([yExtent[0], yExtent[1]])
                .range([height, 0])
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
            .text(`Cell ${this.data.cell_num}: Gene Expression Analysis`);

        this.xLabel = this.svg.append('text')
            .attr('class', 'x-label')
            .attr('text-anchor', 'middle')
            .attr('x', this.width / 2)
            .attr('y', this.height + PLOT_CONFIG.margin.bottom - 5)
            .text('pciSeq: Cell Gene Counts');

        this.yLabel = this.svg.append('text')
            .attr('class', 'y-label')
            .attr('text-anchor', 'middle')
            .attr('transform', 'rotate(-90)')
            .attr('x', -this.height / 2)
            .attr('y', -PLOT_CONFIG.margin.left + 15);

        this.updateLabels();
    }

    updateLabels() {
        this.yLabel.text('scRNAseq: Class Gene Counts');
        
        // this.subtitle.text('Comparing observed counts with expected counts for selected cell type');
    }

    setupClassSelector() {
        // Move dropdown further up to avoid title overlap
        const dropdownGroup = this.svg.append("g")
            .attr("class", "dropdown-group")
            .attr("transform", `translate(0, -85)`);  // Changed from -45 to -65

        // Add background rectangle
        dropdownGroup.append("rect")
            .attr("width", 200)
            .attr("height", 30)
            .attr("fill", "white")
            .attr("opacity", 0.8)
            .attr("rx", 4);

        // Create array of class options with probabilities
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
            
            // Update labels first
            this.updateLabels();
            
            // Then update the plot
            this.updatePlot(false);
        });
    }

    setup() {
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        this.svg = this.createSvg();
        this.scales = this.createScales(width, height);

        this.setupAxes();
        this.setupLabels();
        this.setupDiagonalLine();
        this.setupClassSelector();
    }

    setupDiagonalLine() {
        this.diagonalLine = this.svg.append('line')
            .attr('class', 'diagonal-line')
            .style('stroke', 'gray')
            .style('stroke-width', '1px')
            .style('stroke-dasharray', '4');
        this.updateDiagonalLine();
    }

    setupRegressionLine() {
        this.regressionLine = this.svg.append('line')
            .attr('class', 'regression-line')
            .style('stroke', '#FF9999')
            .style('stroke-width', '2px')
            .style('stroke-dasharray', '3,3')
            .style('opacity', 0);
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
            .style("gap", "5px");

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
            .text("Show regression line");

        this.rSquaredText = div.append("span")
            .style("opacity", 0)
            .style("margin-left", "5px");
    }

    calculateRegression(plotData) {
        const x = plotData.map(d => d.x);
        const y = plotData.map(d => d.y);

        const xMean = d3.mean(x);
        const yMean = d3.mean(y);

        let ssxx = 0, ssyy = 0, ssxy = 0;
        for (let i = 0; i < x.length; i++) {
            ssxx += (x[i] - xMean) * (x[i] - xMean);
            ssyy += (y[i] - yMean) * (y[i] - yMean);
            ssxy += (x[i] - xMean) * (y[i] - yMean);
        }

        const slope = ssxy / ssxx;
        const intercept = yMean - slope * xMean;
        const rSquared = (ssxy * ssxy) / (ssxx * ssyy);

        return { slope, intercept, rSquared };
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
            .attr('cx', d => this.scales.x(d.x))
            .attr('cy', d => this.scales.y(d.y_adj))
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
                this.tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);

                this.tooltip.html(
                    `<strong>${d.name}</strong><br>` +
                    `Observed: ${d.x.toFixed(2)}<br>` +
                    `Expected (raw): ${(d.y).toFixed(2)}<br>` +  // Show raw expected
                    `Expected (adjusted): ${d.y_adj.toFixed(2)}<br>` +  // Show efficiency-adjusted
                    `Efficiency: ${efficiency.toFixed(2)}`
                );

                const tooltipWidth = this.tooltip.node().getBoundingClientRect().width;
                let left = event.pageX;

                if (event.pageX + tooltipWidth/2 > window.innerWidth) {
                    left = window.innerWidth - tooltipWidth - 10;
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

        this.rSquaredText
            .style('opacity', 1)
            .text(`R² = ${rSquared.toFixed(3)}`);
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

        this.diagonalLine
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('x1', this.scales.x(minDomain))
            .attr('y1', this.scales.y(minDomain))
            .attr('x2', this.scales.x(maxDomain))
            .attr('y2', this.scales.y(maxDomain));
    }

    updatePlot(initial = false) {
        // Get index of current class
        const classIndex = this.data.class_names.indexOf(this.currentUserClass);

        const plotData = this.data.gene_names
            .map((name, i) => ({
                name,
                g: +this.data.gene_efficiency[i],
                x: this.data.gene_counts[i],
                y: this.data.scRNAseq_gene_counts[i][classIndex],
                y_adj: this.data.scRNAseq_gene_counts[i][classIndex] * this.data.gene_efficiency[i]
            }))
            .filter(d => this.visibleGenes.has(d.name));

        // Only recalculate scales on initial render
        if (initial) {
            this.scales = this.createScales(this.width, this.height);

            this.xAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisBottom(this.scales.x));

            this.yAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisLeft(this.scales.y));
        } else {
            // For class changes, only update y-axis scale
            const yValues = this.data.scRNAseq_gene_counts.map(row => row[classIndex]);
            const yValues_adj = yValues.map((value, index) => value * this.data.gene_efficiency[index]);
            const yExtent = d3.extent(yValues_adj);
            const yPadding = (yExtent[1] - yExtent[0]) * 0.05;

            this.scales.y.domain([yExtent[0], yExtent[1]]);

            this.yAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(d3.axisLeft(this.scales.y));
        }

        // Update title to show selected class
        this.title.text(`Cell ${this.data.cell_num}: ${this.currentUserClass} - Observed vs Expected Counts`);

        this.updatePoints(plotData);
        this.updateDiagonalLine();
    }

    updateVisibleGenes(visibleGenes) {
        this.visibleGenes = visibleGenes;
        this.updatePlot(false);
    }

    resize() {
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        this.svg
            .attr('width', width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom);

        this.scales.x.range([0, width]);
        this.scales.y.range([height, 0]);

        this.updatePlot(false);
    }
}

