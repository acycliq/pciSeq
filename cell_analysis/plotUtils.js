import { PLOT_CONFIG } from './plotConfig.js';

export function formatProbability(prob) {
    return `${(prob * 100).toFixed(2)}%`;
}

export function createScales(width, height, xData, yData) {
    const xExtent = d3.extent(xData);
    const yExtent = d3.extent(yData);

    // Add a small padding to the domain
    const xPadding = (xExtent[1] - xExtent[0]) * 0.05;
    const yPadding = (yExtent[1] - yExtent[0]) * 0.05;

    return {
        x: d3.scaleLinear()
            .domain([xExtent[0] - xPadding, xExtent[1] + xPadding])
            .range([0, width]),
        y: d3.scaleLinear()
            .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
            .range([height, 0])
    };
}

export function createAxis(scale, orientation) {
    const axis = orientation === 'bottom' ? d3.axisBottom(scale) : d3.axisLeft(scale);
    return axis
        .tickSize(PLOT_CONFIG.axis.tickSize)
        .tickPadding(PLOT_CONFIG.axis.tickPadding);
}
