// Save this as: static/js/main.js

import { ScatterPlot } from './scatterPlot.js';

export function initializePlots(data) {
    // Create main plot instance
    const mainPlot = new ScatterPlot('plot-area', data);

    // Create bottom plot instance
    const bottomPlot = new ScatterPlot('bottom-plot', data);

    // Initialize class selector
    setupClassSelector(data, mainPlot);

    // Add window resize handler
    window.addEventListener('resize', () => {
        mainPlot.resize();
        bottomPlot.resize();
    });

    return { mainPlot, bottomPlot };
}

function setupClassSelector(data, plot) {
    // Get the selector element
    const selector = document.getElementById('class-selector');

    // Prepare options for the dropdown
    const options = prepareSelectorOptions(
        data.class_names,
        data.class_probs,
        data.assigned_class
    );

    // Clear existing options
    selector.innerHTML = '';

    // Add options to selector
    options.forEach(option => {
        const optElement = document.createElement('option');
        optElement.value = option.value;
        optElement.textContent = option.label;
        optElement.selected = option.value === data.user_class;
        selector.appendChild(optElement);
    });

    // Add change event listener
    selector.addEventListener('change', (event) => {
        plot.currentUserClass = event.target.value;
        plot.updatePlot();
    });
}

function prepareSelectorOptions(classes, classProbs, currentAssignedClass) {
    // Filters out the assigned class and formats the dropdown options
    return classes
        .filter(c => c !== currentAssignedClass)
        .map(c => ({
            value: c,
            label: `${c} (${(classProbs[c] * 100).toFixed(2)}%)`,
            probability: classProbs[c]
        }));
}
