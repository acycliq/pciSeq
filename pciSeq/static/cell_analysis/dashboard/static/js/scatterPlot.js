import { PLOT_CONFIG, calculateDimensions } from './plotConfig.js';
import { InterpretationGuide } from './interpretationGuide.js';

export class ScatterPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.currentUserClass = data.user_class;
        this.currentAssignedClass = data.assigned_class;
        this.tooltip = tooltip;
        this.sizeByCount = false;
        this.defaultRadius = 5.5;
        this.visibleGenes = new Set(data.gene_names);

        this.radiusScale = d3.scaleSqrt()
            .domain([0, d3.max(data.gene_counts) || 1])
            .range([2, 8]);

        this.setup();
        this.updatePlot(true);
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
        return geneNames
            .map((gene, index) => ({
                name: gene,
                x: contrData[assignedClass][index],
                y: contrData[userClass][index],
                geneCount: this.data.gene_counts[index] || 0
            }))
            .filter(d => this.visibleGenes.has(d.name));
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
    }

    createSvg() {
        return d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', this.width + PLOT_CONFIG.margin.left + PLOT_CONFIG.margin.right)
            .attr('height', this.height + PLOT_CONFIG.margin.top + PLOT_CONFIG.margin.bottom)
            .append('g')
            .attr('transform', `translate(${PLOT_CONFIG.margin.left},${PLOT_CONFIG.margin.top})`);
    }

    setupAxes() {
        this.xAxis = this.svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${this.height})`)
            .call(this.createAxis(this.scales.x, 'bottom'));

        this.yAxis = this.svg.append('g')
            .attr('class', 'y-axis')
            .call(this.createAxis(this.scales.y, 'left'));
    }

    setupLabels() {
        this.title = this.svg.append("text")
            .attr("class", "plot-title")
            .attr("x", this.width / 2)
            .attr("y", -PLOT_CONFIG.margin.top / 2)
            .style("text-anchor", "middle")
            .text(`Loglikelihood contributions for cell num: ${this.data.cell_num}`);

        this.subtitle = this.svg.append("text")
            .attr("class", "plot-subtitle")
            .attr("x", this.width / 2)
            .attr("y", -PLOT_CONFIG.margin.top / 2 + 20)
            .style("text-anchor", "middle")
            .text(`Assigned class: ${this.currentAssignedClass}`);

        this.xLabel = this.svg.append("text")
            .attr("class", "axis-label x-axis-label")
            .attr("x", this.width / 2)
            .attr("y", this.height + PLOT_CONFIG.margin.bottom - 10)
            .style("text-anchor", "middle")
            .text(this.getAxisLabel(this.currentAssignedClass, this.data.class_probs[this.currentAssignedClass]));

        this.yLabel = this.svg.append("text")
            .attr("class", "axis-label y-axis-label")
            .attr("transform", "rotate(-90)")
            .attr("x", -this.height / 2)
            .attr("y", -PLOT_CONFIG.margin.left + 40)
            .style("text-anchor", "middle")
            .text(this.getAxisLabel(this.currentUserClass, this.data.class_probs[this.currentUserClass]));
    }

    setupDiagonalLine() {
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

    setupSizeControl() {
        const checkboxGroup = this.svg.append("g")
            .attr("class", "checkbox-group")
            .attr("transform", `translate(${this.width - 200}, 20)`);

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
                this.updatePlot(false);
            });
    }

    getAxisLabel(className, probability) {
        return `${className} (Prob: ${(probability * 100).toFixed(2)}%)`;
    }

    updateVisibleGenes(visibleGenes) {
        this.visibleGenes = new Set(visibleGenes);
        this.updatePlot(false);
    }

    updatePlot(updateScales = true) {
        const plotData = this.preparePlotData(
            this.data.gene_names,
            this.data.contr,
            this.currentAssignedClass,
            this.currentUserClass
        );

        if (updateScales) {
            this.scales = this.createScales(
                this.width,
                this.height,
                this.data.contr[this.currentAssignedClass],
                this.data.contr[this.currentUserClass]
            );

            this.xAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(this.createAxis(this.scales.x, 'bottom'));

            this.yAxis.transition()
                .duration(PLOT_CONFIG.animation.duration)
                .call(this.createAxis(this.scales.y, 'left'));
        }

        this.updateLabels();
        this.updatePoints(plotData);
        this.updateDiagonalLine();
        this.guide.update(this.currentUserClass, this.currentAssignedClass);
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
            .attr('cy', d => this.scales.y(d.y))
            .attr('r', d => this.sizeByCount ?
                this.radiusScale(d.geneCount) :
                this.defaultRadius);

        this.svg.selectAll('circle')
            .on('mouseenter', (event, d) => {
                d3.select(event.target)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeIn)
                    .attr('r', d => this.sizeByCount ?
                        this.radiusScale(d.geneCount) * 1.6 :
                        this.defaultRadius * 1.6)
                    .attr('fill', PLOT_CONFIG.point.hoverColor);
            })
            .on('mouseleave', (event, d) => {
                d3.select(event.target)
                    .transition()
                    .duration(PLOT_CONFIG.animation.tooltip.fadeOut)
                    .attr('r', d => this.sizeByCount ?
                        this.radiusScale(d.geneCount) :
                        this.defaultRadius)
                    .attr('fill', PLOT_CONFIG.point.color);
            })
            .on('mouseover', (event, d) => {
                const vw = Math.max(document.documentElement.clientWidth || 0, window.innerWidth || 0);

                this.tooltip.transition()
                    .duration(200)
                    .style("opacity", .9);

                this.tooltip.html(
                    `<strong>${d.name}</strong><br>` +
                    `X: ${d.x.toFixed(3)}<br>` +
                    `Y: ${d.y.toFixed(3)}<br>` +
                    `Gene Count: ${d.geneCount.toFixed(2)}`
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

    updateLabels() {
        this.yLabel.text(this.getAxisLabel(
            this.currentUserClass,
            this.data.class_probs[this.currentUserClass]
        ));

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