/**
 * Setup Module
 * Handles configuration form, file selection, and job launching
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};

    // Setup state
    const setupState = {
        isRunning: false
    };

    // Initialize setup form handlers
    function initializeSetupForm() {
        // File input handlers
        setupFileHandlers();

        // Form submission
        const form = document.getElementById('setup-form');
        if (form) {
            form.addEventListener('submit', handleFormSubmit);
        }
    }

    // Setup file input handlers (browser + paste)
    function setupFileHandlers() {
        // Spots file
        const spotsFile = document.getElementById('spots-file');
        const spotsPath = document.getElementById('spots-path');
        if (spotsFile && spotsPath) {
            spotsFile.addEventListener('change', (e) => {
                if (e.target.files[0]) {
                    spotsPath.value = e.target.files[0].path || e.target.files[0].name;
                }
            });
        }

        // scRNAseq file
        const scrnaFile = document.getElementById('scrna-file');
        const scrnaPath = document.getElementById('scrna-path');
        if (scrnaFile && scrnaPath) {
            scrnaFile.addEventListener('change', (e) => {
                if (e.target.files[0]) {
                    scrnaPath.value = e.target.files[0].path || e.target.files[0].name;
                }
            });
        }

        // Cell masks file
        const cooFile = document.getElementById('coo-file');
        const cooPath = document.getElementById('coo-path');
        if (cooFile && cooPath) {
            cooFile.addEventListener('change', (e) => {
                if (e.target.files[0]) {
                    cooPath.value = e.target.files[0].path || e.target.files[0].name;
                }
            });
        }
    }

    // Validate form inputs
    function validateForm() {
        let isValid = true;
        const errors = {};

        // Required fields
        const spotsPath = document.getElementById('spots-path').value.trim();
        const scrnaPath = document.getElementById('scrna-path').value.trim();
        const cooPath = document.getElementById('coo-path').value.trim();

        if (!spotsPath) {
            errors.spots = 'Spots CSV path is required';
            isValid = false;
        }

        if (!scrnaPath) {
            errors.scrna = 'scRNAseq CSV path is required';
            isValid = false;
        }

        if (!cooPath) {
            errors.coo = 'Cell masks path is required';
            isValid = false;
        }

        // Validate voxel size format (x,y,z)
        const voxelSize = document.getElementById('param-voxelsize').value.trim();
        const voxelParts = voxelSize.split(',').map(v => v.trim());
        if (voxelParts.length !== 3 || voxelParts.some(v => isNaN(parseFloat(v)))) {
            errors.voxelsize = 'Voxel size must be in format: x,y,z (e.g., 1,1,1)';
            isValid = false;
        }

        // Display errors
        displayErrors(errors);

        return isValid;
    }

    // Display validation errors
    function displayErrors(errors) {
        // Clear all errors first
        ['spots', 'scrna', 'coo'].forEach(field => {
            const errorEl = document.getElementById(`${field}-error`);
            if (errorEl) errorEl.textContent = '';
        });

        // Show new errors
        Object.entries(errors).forEach(([field, message]) => {
            const errorEl = document.getElementById(`${field}-error`);
            if (errorEl) {
                errorEl.textContent = message;
            }
        });
    }

    // Handle form submission
    async function handleFormSubmit(e) {
        e.preventDefault();

        if (setupState.isRunning) {
            return; // Prevent double submission
        }

        // Validate form
        if (!validateForm()) {
            return;
        }

        // Collect configuration
        const config = collectConfiguration();

        // Disable start button
        setupState.isRunning = true;
        const startBtn = document.getElementById('start-btn');
        if (startBtn) {
            startBtn.disabled = true;
            startBtn.textContent = 'Starting...';
        }

        try {
            // Send configuration to backend
            const response = await fetch('/api/start_job', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify(config)
            });

            // Check if response is JSON
            const contentType = response.headers.get('content-type');
            if (!contentType || !contentType.includes('application/json')) {
                const text = await response.text();
                console.error('Server returned non-JSON response:', text);
                throw new Error(`Server error: ${response.status} ${response.statusText}. Check server console for details.`);
            }

            const result = await response.json();

            if (response.ok) {
                console.log('Job started successfully:', result);
                // Switch to visualization view
                switchToVisualizationView();
            } else {
                throw new Error(result.error || 'Failed to start job');
            }
        } catch (error) {
            console.error('Error starting job:', error);
            alert(`Error starting analysis: ${error.message}`);

            // Re-enable start button
            setupState.isRunning = false;
            if (startBtn) {
                startBtn.disabled = false;
                startBtn.textContent = 'Start Analysis';
            }
        }
    }

    // Collect configuration from form
    function collectConfiguration() {
        // Parse voxel size
        const voxelSizeStr = document.getElementById('param-voxelsize').value.trim();
        const voxelSize = voxelSizeStr.split(',').map(v => parseFloat(v.trim()));

        return {
            // Data inputs
            spots_path: document.getElementById('spots-path').value.trim(),
            scrna_path: document.getElementById('scrna-path').value.trim(),
            coo_path: document.getElementById('coo-path').value.trim(),

            // Parameters
            Inefficiency: parseFloat(document.getElementById('param-inefficiency').value),
            nNeighbors: parseInt(document.getElementById('param-nneighbors').value),
            CellCallTolerance: parseFloat(document.getElementById('param-tolerance').value),
            rSpot: parseFloat(document.getElementById('param-rspot').value),
            max_iter: parseInt(document.getElementById('param-maxiter').value),
            MisreadDensity: parseFloat(document.getElementById('param-misread').value),
            voxel_size: voxelSize,
            output_path: document.getElementById('param-output').value.trim(),

            // Options
            save_data: document.getElementById('opt-save').checked,
            remove_flat_cells: document.getElementById('opt-remove-flat').checked,
            launch_diagnostics: document.getElementById('opt-diagnostics').checked,
            realtime_viewer: true, // Always true since we're in the viewer
            realtime_viewer_port: 5001
        };
    }

    // Switch from setup view to visualization view
    function switchToVisualizationView() {
        const setupView = document.getElementById('setup-view');
        const visualizationView = document.getElementById('container');

        if (setupView) {
            setupView.classList.add('view-hidden');
        }

        if (visualizationView) {
            visualizationView.classList.remove('view-hidden');
        }

        // Initialize visualization components (if not already initialized)
        // The socket handlers will automatically connect and start receiving data
        console.log('Switched to visualization view');
    }

    // Export functions
    window.pciSeq.setup = {
        initializeSetupForm: initializeSetupForm,
        switchToVisualizationView: switchToVisualizationView
    };

})();