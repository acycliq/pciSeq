// Save this as: static/js/index.js

import { initializePlots } from './main.js';
import { DataService } from './api/dataService.js';

// Initialize the application
async function init() {
    try {
        // Show loading state
        document.getElementById('plot-area').innerHTML = 'Loading...';
        document.getElementById('bottom-plot').innerHTML = 'Loading...';

        // Fetch data
        const urlParams = new URLSearchParams(window.location.search);
        const cellId = urlParams.get('cell_id') || '1'; // Default to cell 1 if not specified
        const data = await DataService.fetchCellData(cellId);

        // Clear loading state
        document.getElementById('plot-area').innerHTML = '';
        document.getElementById('bottom-plot').innerHTML = '';

        // Initialize plots
        const plots = initializePlots(data);

        // Make plots available globally for debugging
        window.plots = plots;

    } catch (error) {
        // Handle errors
        console.error('Initialization failed:', error);
        document.getElementById('plot-area').innerHTML =
            `<div class="error">Failed to load plot: ${error.message}</div>`;
        document.getElementById('bottom-plot').innerHTML =
            `<div class="error">Failed to load plot: ${error.message}</div>`;
    }
}

// Start the application when DOM is ready
document.addEventListener('DOMContentLoaded', init);
