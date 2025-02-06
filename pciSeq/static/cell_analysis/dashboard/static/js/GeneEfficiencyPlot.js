import { PLOT_CONFIG } from './plotConfig.js';

export class GeneEfficiencyPlot {
    constructor(containerId, data, tooltip) {
        this.containerId = containerId;
        this.data = data;
        this.tooltip = tooltip;
        this.margin = PLOT_CONFIG.margin;
        this.visibleGenes = new Set(data.gene_names);
        
        // Create SVG
        this.svg = d3.select(`#${this.containerId}`)
            .append("svg")
            .style("width", "100%")
            .style("height", "100%");
            
        // Add group for plot
        this.g = this.svg.append("g");
        
        // Create axes groups
        this.xAxis = this.g.append("g")
            .attr("class", "x-axis")
            .style("font-size", PLOT_CONFIG.axis.fontSize);
            
        this.yAxis = this.g.append("g")
            .attr("class", "y-axis")
            .style("font-size", PLOT_CONFIG.axis.fontSize);
        
        // Add title
        this.svg.append("text")
            .attr("class", "plot-title")
            .attr("text-anchor", "middle")
            .style("font-size", "18px")
            .style("font-weight", "600")
            .text("Gene Inefficiency");
            
        this.resizeHandler = () => this.resize();
        window.addEventListener('resize', this.resizeHandler);
        
        this.resize();
        this.updatePlot(true);
    }

    resize() {
        const bbox = document.getElementById(this.containerId).getBoundingClientRect();
        this.width = bbox.width - this.margin.left - this.margin.right;
        this.height = bbox.height - this.margin.top - this.margin.bottom;

        this.svg.attr("viewBox", `0 0 ${bbox.width} ${bbox.height}`);
        this.g.attr("transform", `translate(${this.margin.left},${this.margin.top})`);
        
        this.svg.select(".plot-title")
            .attr("x", bbox.width / 2)
            .attr("y", this.margin.top / 2);

        this.updatePlot(true);
    }

    updateVisibleGenes(visibleGenes) {
        this.visibleGenes = new Set(visibleGenes);
        this.updatePlot(false);
    }

    updatePlot(initial = false) {
        if (!this.data || !this.width) return;

        const plotData = this.data.gene_names
            .map((gene, i) => ({
                gene: gene,
                efficiency: this.data.gene_efficiency[i]
            }))
            .filter(d => this.visibleGenes.has(d.gene))
            .sort((a, b) => a.gene.localeCompare(b.gene));

        // Calculate max efficiency and add 5% padding
        const maxEfficiency = d3.max(plotData, d => d.efficiency);
        const padding = maxEfficiency * 0.05;

        const xScale = d3.scaleBand()
            .domain(plotData.map(d => d.gene))
            .range([0, this.width])
            .padding(0.1);

        const yScale = d3.scaleLinear()
            .domain([0, maxEfficiency + padding])
            .range([this.height, 0]);

        // Update axes with transitions
        this.xAxis
            .attr("transform", `translate(0,${this.height})`)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .call(d3.axisBottom(xScale)
                .tickFormat((d, i) => {
                    return i % 10 === 0 ? d : '';
                }));
            
        // Apply text rotation after transition
        this.xAxis.selectAll("text")
            .style("text-anchor", "end")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", "rotate(-45)");

        this.yAxis
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .call(d3.axisLeft(yScale).ticks(5));

        // Update points
        const points = this.g.selectAll("circle")
            .data(plotData, d => d.gene);

        points.exit()
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr("r", 0)
            .remove();

        const pointsEnter = points.enter()
            .append("circle")
            .attr("r", 0)
            .attr("cx", d => xScale(d.gene) + xScale.bandwidth()/2)  // Center in band
            .attr("cy", d => yScale(d.efficiency));

        points.merge(pointsEnter)
            .transition()
            .duration(PLOT_CONFIG.animation.duration)
            .attr("cx", d => xScale(d.gene) + xScale.bandwidth()/2)
            .attr("cy", d => yScale(d.efficiency))
            .attr("r", PLOT_CONFIG.point.radius)
            .attr("fill", PLOT_CONFIG.point.color)
            .attr("stroke", "white")
            .attr("stroke-width", "0.5");

        // Hover effects remain the same
        this.g.selectAll("circle")
            .on("mouseover", (event, d) => {
                d3.select(event.currentTarget)
                    .transition()
                    .duration(100)
                    .attr("r", PLOT_CONFIG.point.radius * 1.5)
                    .style("fill", PLOT_CONFIG.point.hoverColor);
                    
                this.tooltip
                    .style("opacity", 0.9)
                    .html(`<strong>${d.gene}</strong><br>Inefficiency: ${d.efficiency.toFixed(3)}`)
                    .style("left", `${event.pageX + 10}px`)
                    .style("top", `${event.pageY - 28}px`);
            })
            .on("mouseout", (event) => {
                d3.select(event.currentTarget)
                    .transition()
                    .duration(100)
                    .attr("r", PLOT_CONFIG.point.radius)
                    .style("fill", PLOT_CONFIG.point.color);
                    
                this.tooltip
                    .style("opacity", 0);
            });
    }

    destroy() {
        window.removeEventListener('resize', this.resizeHandler);
    }
}