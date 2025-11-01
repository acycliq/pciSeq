/**
 * Convergence Chart Module
 * D3.js-based visualization of algorithm convergence over iterations
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Chart state
    let chartSvg = null;
    let chartScales = { x: null, y: null };
    let chartTooltip = null;

    // Render convergence chart
    function renderConvergenceChart() {
        const container = document.getElementById('convergence-svg-container');
        if (!container) return;

        const containerWidth = container.clientWidth;
        const containerHeight = container.clientHeight;

        // Margins for axes and labels
        const margin = { top: 15, right: 15, bottom: 30, left: 45 };
        const width = containerWidth - margin.left - margin.right;
        const height = containerHeight - margin.top - margin.bottom;

        // Initialize SVG on first call
        if (!chartSvg) {
            // Create tooltip
            chartTooltip = d3.select('body')
                .append('div')
                .attr('class', 'chart-tooltip');

            // Create SVG
            chartSvg = d3.select(container)
                .append('svg')
                .attr('width', '100%')
                .attr('height', '100%')
                .append('g')
                .attr('transform', `translate(${margin.left},${margin.top})`);

            // Add grid group (drawn first, behind everything)
            chartSvg.append('g').attr('class', 'grid');

            // Add threshold line group
            chartSvg.append('g').attr('class', 'threshold-line');

            // Add line path group
            chartSvg.append('g').attr('class', 'line-group');

            // Add dots group
            chartSvg.append('g').attr('class', 'dots');

            // Add axes groups
            chartSvg.append('g').attr('class', 'x-axis');
            chartSvg.append('g').attr('class', 'y-axis');

            // Add axis labels
            chartSvg.append('text')
                .attr('class', 'x-axis-label')
                .attr('text-anchor', 'middle')
                .attr('x', width / 2)
                .attr('y', height + 28)
                .style('fill', '#e5e5e5')
                .style('font-size', '11px')
                .text('Iteration');

            chartSvg.append('text')
                .attr('class', 'y-axis-label')
                .attr('text-anchor', 'middle')
                .attr('transform', `translate(-35, ${height / 2}) rotate(-90)`)
                .style('fill', '#e5e5e5')
                .style('font-size', '11px')
                .text('Delta (Convergence)');
        }

        // No data yet
        if (state.deltaHistory.length === 0) {
            chartSvg.selectAll('.waiting-text').remove();
            chartSvg.append('text')
                .attr('class', 'waiting-text')
                .attr('x', width / 2)
                .attr('y', height / 2)
                .attr('text-anchor', 'middle')
                .style('fill', '#9aa0a6')
                .style('font-size', '11px')
                .text('Waiting for data...');
            return;
        } else {
            chartSvg.selectAll('.waiting-text').remove();
        }

        // Update scales
        const xScale = d3.scaleLinear()
            .domain([0, state.chartXMax])
            .range([0, width]);

        const yScale = d3.scaleLinear()
            .domain([0, 1.0])
            .range([height, 0]);

        chartScales.x = xScale;
        chartScales.y = yScale;

        // Update grid
        const grid = chartSvg.select('.grid');
        grid.selectAll('*').remove();

        // Horizontal grid lines
        for (let y = 0; y <= 1.0; y += 0.2) {
            grid.append('line')
                .attr('x1', 0)
                .attr('x2', width)
                .attr('y1', yScale(y))
                .attr('y2', yScale(y))
                .style('stroke', '#2a2a2a')
                .style('stroke-width', 1);
        }

        // Vertical grid lines
        const xStep = state.chartXMax <= 100 ? 25 : 50;
        for (let x = 0; x <= state.chartXMax; x += xStep) {
            grid.append('line')
                .attr('x1', xScale(x))
                .attr('x2', xScale(x))
                .attr('y1', 0)
                .attr('y2', height)
                .style('stroke', '#2a2a2a')
                .style('stroke-width', 1);
        }

        // Update threshold line
        const threshold = chartSvg.select('.threshold-line');
        threshold.selectAll('*').remove();

        threshold.append('line')
            .attr('x1', 0)
            .attr('x2', width)
            .attr('y1', yScale(state.cellCallTolerance))
            .attr('y2', yScale(state.cellCallTolerance))
            .style('stroke', '#FFD54F')
            .style('stroke-width', 1.5)
            .style('stroke-dasharray', '4,4');

        threshold.append('text')
            .attr('x', width - 2)
            .attr('y', yScale(state.cellCallTolerance) - 3)
            .attr('text-anchor', 'end')
            .style('fill', '#FFD54F')
            .style('font-size', '10px')
            .text(`Tolerance (${state.cellCallTolerance.toFixed(2)})`);

        // Update line path
        const line = d3.line()
            .x(d => xScale(d.iteration))
            .y(d => yScale(d.delta))
            .curve(d3.curveMonotoneX); // Smooth curve

        const lineGroup = chartSvg.select('.line-group');
        const path = lineGroup.selectAll('.convergence-line').data([state.deltaHistory]);

        path.enter()
            .append('path')
            .attr('class', 'convergence-line')
            .style('fill', 'none')
            .style('stroke', '#22c55e')
            .style('stroke-width', 2)
            .merge(path)
            .transition()
            .duration(300)
            .ease(d3.easeLinear)
            .attr('d', line);

        // Update dots
        const dots = chartSvg.select('.dots')
            .selectAll('.dot')
            .data(state.deltaHistory, d => d.iteration);

        // Enter new dots with animation
        dots.enter()
            .append('circle')
            .attr('class', 'dot')
            .attr('cx', d => xScale(d.iteration))
            .attr('cy', d => yScale(d.delta))
            .attr('r', 0)
            .style('fill', '#22c55e')
            .style('cursor', 'pointer')
            .on('mouseover', function(event, d) {
                d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('r', 5);

                chartTooltip
                    .classed('visible', true)
                    .html(`<strong>Iteration ${d.iteration}</strong><br/>Delta: ${d.delta.toFixed(4)}`)
                    .style('left', (event.pageX + 10) + 'px')
                    .style('top', (event.pageY - 10) + 'px');
            })
            .on('mouseout', function() {
                d3.select(this)
                    .transition()
                    .duration(100)
                    .attr('r', 3);

                chartTooltip.classed('visible', false);
            })
            .transition()
            .duration(300)
            .attr('r', 3);

        // Update existing dots
        dots.transition()
            .duration(300)
            .attr('cx', d => xScale(d.iteration))
            .attr('cy', d => yScale(d.delta));

        // Remove old dots
        dots.exit().remove();

        // Update axes
        const xAxis = d3.axisBottom(xScale)
            .ticks(state.chartXMax <= 100 ? 5 : 10)
            .tickFormat(d3.format('d'));

        const yAxis = d3.axisLeft(yScale)
            .ticks(5)
            .tickFormat(d3.format('.1f'));

        chartSvg.select('.x-axis')
            .attr('transform', `translate(0, ${height})`)
            .call(xAxis)
            .style('color', '#9aa0a6')
            .style('font-size', '10px');

        chartSvg.select('.y-axis')
            .call(yAxis)
            .style('color', '#9aa0a6')
            .style('font-size', '10px');

        // Style axis lines and ticks
        chartSvg.selectAll('.x-axis line, .y-axis line, .x-axis path, .y-axis path')
            .style('stroke', '#9aa0a6');
    }

    // Export functions
    window.pciSeq.chart = {
        renderConvergenceChart: renderConvergenceChart
    };

})();