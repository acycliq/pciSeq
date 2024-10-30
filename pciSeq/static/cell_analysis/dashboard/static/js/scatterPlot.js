// Save this as: static/js/components/ScatterPlot.js

import { PLOT_CONFIG, calculateDimensions } from './plotConfig.js';
import { InterpretationGuide } from './interpretationGuide.js';

export class ScatterPlot {
    constructor(containerId, data, tooltip) {
        // Initialize the plot with data and container ID
        this.containerId = containerId;
        this.data = data;
        this.currentUserClass = data.user_class;
        this.currentAssignedClass = data.assigned_class;
        this.tooltip = tooltip
        
        this.setup();
        this.updatePlot();
    }

    createScales(width, height, xData, yData) {
        return {
            x: d3.scaleLinear()
                .domain(d3.extent(xData))
                .range([0, width])
                .nice(),
            y: d3.scaleLinear()
                .domain(d3.extent(yData))
                .range([height, 0])
                .nice()
        };
    }

    createAxis(scale, orientation) {
        const axis = orientation === 'bottom' ? d3.axisBottom : d3.axisLeft;
        return axis(scale)
            .ticks(5)
            .tickSize(-5)
            .tickPadding(5);
    }

    preparePlotData(geneNames, contrData, assignedClass, userClass) {
        // Maps gene names to their corresponding x,y coordinates
        return geneNames.map((gene, index) => ({
            name: gene,
            x: contrData[assignedClass][index],
            y: contrData[userClass][index]
        }));
    }


    setup() {
        // Set up the basic plot structure
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
        
        // this.tooltip = createTooltip();
        this.guide = new InterpretationGuide(this.svg, width, height);
        
        this.setupAxes();
        this.setupLabels();
        this.setupDiagonalLine();
    }

    createSvg() {
        // Create the SVG container
        return d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', this.width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', this.height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom)
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top})`);
    }

    setupAxes() {
        // Create and append x-axis
        this.xAxis = this.svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${this.height})`)
            .call(this.createAxis(this.scales.x, 'bottom'));

        // Create and append y-axis
        this.yAxis = this.svg.append('g')
            .attr('class', 'y-axis')
            .call(this.createAxis(this.scales.y, 'left'));
    }

    setupLabels() {
        // Add plot title
        this.title = this.svg.append("text")
            .attr("class", "plot-title")
            .attr("x", this.width / 2)
            .attr("y", -PLOT_CONFIG.margin.top / 2)
            .style("text-anchor", "middle")
            .text(`Loglikelihood contributions for cell num: ${this.data.cell_num}`);

        // Add plot subtitle
        this.subtitle = this.svg.append("text")
            .attr("class", "plot-subtitle")
            .attr("x", this.width / 2)
            .attr("y", -PLOT_CONFIG.margin.top / 2 + 20)
            .style("text-anchor", "middle")
            .text(`Assigned class: ${this.currentAssignedClass}`);

        // Add x-axis label
        this.xLabel = this.svg.append("text")
            .attr("class", "axis-label x-axis-label")
            .attr("x", this.width / 2)
            .attr("y", this.height + PLOT_CONFIG.margin.bottom - 10)
            .style("text-anchor", "middle")
            .text(this.getAxisLabel(this.currentAssignedClass, this.data.class_probs[this.currentAssignedClass]));

        // Add y-axis label
        this.yLabel = this.svg.append("text")
            .attr("class", "axis-label y-axis-label")
            .attr("transform", "rotate(-90)")
            .attr("x", -this.height / 2)
            .attr("y", -PLOT_CONFIG.margin.left + 40)
            .style("text-anchor", "middle")
            .text(this.getAxisLabel(this.currentUserClass, this.data.class_probs[this.currentUserClass]));
    }

    setupDiagonalLine() {
        // Add diagonal reference line (y=x)
        const minDomain = Math.min(
            this.scales.x.domain()[0],
            this.scales.y.domain()[0]
        );
        const maxDomain = Math.max(
            this.scales.x.domain()[1],
            this.scales.y.domain()[1]
        );

        this.diagonalLine = this.svg.append('line')
            .attr('class', 'diagonal-line')
            .attr('x1', this.scales.x(minDomain))
            .attr('y1', this.scales.y(minDomain))
            .attr('x2', this.scales.x(maxDomain))
            .attr('y2', this.scales.y(maxDomain))
            .attr('stroke', PLOT_CONFIG.diagonalLine.color)
            .attr('stroke-width', PLOT_CONFIG.diagonalLine.width)
            .attr('stroke-dasharray', PLOT_CONFIG.diagonalLine.dashArray);
    }

    getAxisLabel(className, probability) {
        return `${className} (Prob: ${(probability * 100).toFixed(2)}%)`;
    }

    updatePlot() {
        // Prepare data for plotting
        const plotData = this.preparePlotData(
            this.data.gene_names,
            this.data.contr,
            this.currentAssignedClass,
            this.currentUserClass
        );

        // Update scales
        this.scales.x.domain(d3.extent(plotData, d => d.x));
        this.scales.y.domain(d3.extent(plotData, d => d.y));

        // Update axes with animation
        this.xAxis.transition()
            .duration(PLOT_CONFIG.animation.duration)
            .call(this.createAxis(this.scales.x, 'bottom'));

        this.yAxis.transition()
            .duration(PLOT_CONFIG.animation.duration)
            .call(this.createAxis(this.scales.y, 'left'));

        // Update labels
        this.updateLabels();

        // Update points
        this.updatePoints(plotData);

        // Update diagonal line
        this.updateDiagonalLine();

        // Update interpretation guide
        this.guide.update(this.currentUserClass, this.currentAssignedClass);
    }

    updatePoints(plotData) {
        // Select all points and bind data
        const dots = this.svg.selectAll('circle')
            .data(plotData);

        // Handle enter selection
        const dotsEnter = dots.enter()
            .append('circle')
            .attr('r', PLOT_CONFIG.point.radius)
            .attr('fill', PLOT_CONFIG.point.color);

        // Handle update selection
        dots.merge(dotsEnter)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr('cx', d => this.scales.x(d.x))
            .attr('cy', d => this.scales.y(d.y))
            .style("stroke", "white")  // add white outline
            .style("stroke-width", 1.5)  // adjust thickness as needed
            .style("stroke-opacity", 0.8);  // slightly transparent outline

        // Handle exit selection
        dots.exit().remove();

        // Add tooltip interactions

        // Add hover events (on the original selection, not the transition)
        this.svg.selectAll('circle')
            .on('mouseenter', function(event, d) {  // Note: need event parameter in D3v6+
                console.log('Mouse enter')
                d3.select(this)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr('r', 1.6 * PLOT_CONFIG.point.radius)
                    .attr('fill', PLOT_CONFIG.point.color);
            })
            .on('mouseleave', function(event, d) {  // Note: need event parameter in D3v6+
                d3.select(this)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr('r', PLOT_CONFIG.point.radius)
                    .attr('fill', PLOT_CONFIG.point.color);
            })
            .on('mouseover', (event, d) => {
                const tooltip = this.tooltip;

                // Get viewport width
                const vw = Math.max(document.documentElement.clientWidth || 0, window.innerWidth || 0);

                // Get tooltip width after setting content
                tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);

                // More compact HTML with bold labels
                tooltip.html(
                    `<strong>${d.name}</strong><br>` +
                    `X: ${d.x.toFixed(3)}<br>` +
                    `Y: ${d.y.toFixed(3)}`
                );

                // Get tooltip dimensions
                const tooltipWidth = tooltip.node().getBoundingClientRect().width;

                // Calculate position
                let left = event.pageX;

                // Adjust if too close to right edge
                if (event.pageX + tooltipWidth/2 > vw) {
                    left = vw - tooltipWidth - 10; // 10px padding from edge
                }
                // Adjust if too close to left edge
                else if (event.pageX - tooltipWidth/2 < 0) {
                    left = tooltipWidth/2 + 10; // 10px padding from edge
                }

                tooltip
                    .style("left", `${left}px`)
                    .style("top", `${event.pageY - 28}px`);
            })
            .on("mouseout", () => {
                this.tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            });
    }

    updateLabels() {
        // Update axis labels with probabilities
        this.yLabel.text(this.getAxisLabel(
            this.currentUserClass,
            this.data.class_probs[this.currentUserClass]
        ));

        this.xLabel.text(this.getAxisLabel(
            this.currentAssignedClass,
            this.data.class_probs[this.currentAssignedClass]
        ));

        // Update subtitle
        this.subtitle.text(
            `Assigned class: ${this.currentAssignedClass} vs Selected class: ${this.currentUserClass}`
        );
    }

    updateDiagonalLine() {
        // Update diagonal line position based on new scales
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

    resize() {
        // Update dimensions
        const { width, height } = calculateDimensions();
        this.width = width;
        this.height = height;

        // Update SVG size
        this.svg.attr('width', width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom);

        // Update scales
        this.scales.x.range([0, width]);
        this.scales.y.range([height, 0]);

        // Trigger complete plot update
        this.updatePlot();
    }
}

