/**
 * Socket.IO Event Handlers Module
 * Manages real-time communication with the pciSeq server
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Socket.IO connection
    const socket = io();

    // Connection events
    socket.on('connect', () => {
        console.log('Connected to pciSeq server');
        window.pciSeq.updateConnectionStatus(true);
    });

    socket.on('disconnect', (reason) => {
        console.log('Disconnected:', reason);
        window.pciSeq.updateConnectionStatus(false);
    });

    // Backward-compatible single-shot update (geometry + classes together)
    socket.on('iteration_update', (data) => {
        console.log('=== RECEIVED iteration_update (compat) ===');
        try {
            if (data.centroids_x && data.centroids_y && data.radii) {
                state.geom.ready = true;
                state.geom.centroids_x = new Float32Array(data.centroids_x);
                state.geom.centroids_y = new Float32Array(data.centroids_y);
                state.geom.radii = new Float32Array(data.radii);
            }
            applyClassesUpdate(data.iteration, data.delta, data.num_cells, data.cell_classes, data.prob);
        } catch (err) {
            console.error('ERROR compat iteration_update:', err);
        }
    });

    // Geometry init (sent once)
    socket.on('geometry_init_begin', (meta) => {
        console.log('=== BEGIN geometry_init ===', meta);
        state.geom.stream = {
            total: meta.num_cells,
            received: 0,
            cell_ids: new Int32Array(meta.num_cells),
            centroids_x: new Float32Array(meta.num_cells),
            centroids_y: new Float32Array(meta.num_cells),
            centroids_z: new Float32Array(meta.num_cells),
            radii: new Float32Array(meta.num_cells)
        };

        // Store meta
        if (meta.mcr !== undefined && meta.mcr !== null) {
            state.geom.mcr = Number(meta.mcr);
        }
        state.geom.is3D = !!meta.is_3d;
        if (meta.voxel_size && Array.isArray(meta.voxel_size) && meta.voxel_size.length === 3) {
            state.geom.voxelSize = meta.voxel_size.map(Number);
        }
        if (meta.img_dim && typeof meta.img_dim === 'object') {
            state.geom.imgDim = meta.img_dim;
        }
        // Store cell_call_tolerance for convergence chart
        if (meta.cell_call_tolerance !== undefined && meta.cell_call_tolerance !== null) {
            state.cellCallTolerance = Number(meta.cell_call_tolerance);
            console.log(`Cell call tolerance set to: ${state.cellCallTolerance}`);
        }
        // Store version (if provided) and update UI
        if (meta.version) {
            state.version = meta.version;
            const versionEl = document.getElementById('version-display');
            if (versionEl) versionEl.textContent = String(state.version);
        }

        // Enable/disable Plane toggle based on 3D meta
        const planeToggleEl = document.getElementById('plane-toggle');
        if (planeToggleEl) {
            planeToggleEl.disabled = !state.geom.is3D;
        }

        // Reset change tracking at the start of a new run/session
        state.previousProb = null;
        state.changedCells.clear();
        window.pciSeq.updateChangesCount();
        const prevNotice = document.getElementById('user-notice');
        if (prevNotice) prevNotice.style.display = 'none';

        // Store class names mapping (index to name)
        if (meta.class_names) {
            state.cellClassNames = {};
            meta.class_names.forEach((name, idx) => {
                state.cellClassNames[idx] = name;
            });
            console.log('Loaded class names:', state.cellClassNames);

            // Apply pending color scheme if one was loaded before class names arrived
            if (state.pendingColorScheme) {
                console.log('Applying pending color scheme...');
                window.pciSeq.colors.applyColorScheme(state.pendingColorScheme);
                state.pendingColorScheme = null;
            }
        }
    });

    socket.on('geometry_init_chunk', (chunk) => {
        if (!state.geom.stream) return;
        const {start, end} = chunk;
        state.geom.stream.cell_ids.set(chunk.cell_ids, start);
        state.geom.stream.centroids_x.set(chunk.centroids_x, start);
        state.geom.stream.centroids_y.set(chunk.centroids_y, start);
        if (chunk.centroids_z) {
            state.geom.stream.centroids_z.set(chunk.centroids_z, start);
        }
        state.geom.stream.radii.set(chunk.radii, start);
        state.geom.stream.received = end;
    });

    socket.on('geometry_init_end', () => {
        console.log('=== END geometry_init ===');
        const gs = state.geom.stream;
        state.geom.cell_ids = gs.cell_ids;
        state.geom.centroids_x = gs.centroids_x;
        state.geom.centroids_y = gs.centroids_y;
        state.geom.centroids_z = gs.centroids_z;
        state.geom.radii = gs.radii;
        state.geom.ready = true;
        state.geom.stream = null;

        console.log(`Geometry arrays: cell_ids=${state.geom.cell_ids ? state.geom.cell_ids.length : 'undefined'}, centroids_x=${state.geom.centroids_x.length}, radii=${state.geom.radii.length}`);

        // Compute plane indices if 3D
        if (state.geom.is3D && state.geom.voxelSize && state.geom.centroids_z) {
            try {
                const dx = Number(state.geom.voxelSize[0]);
                const dz = Number(state.geom.voxelSize[2]);
                const Sz = dz / dx; // isotropic z = original_z * (dz/dx)
                const N = state.geom.centroids_z.length;
                const planeId = new Uint16Array(N);
                let minP = Infinity, maxP = -Infinity;
                for (let i = 0; i < N; i++) {
                    const p = Math.floor(state.geom.centroids_z[i] / Sz);
                    planeId[i] = p;
                    if (p < minP) minP = p;
                    if (p > maxP) maxP = p;
                }
                state.geom.planeId = planeId;
                state.planeRange = { min: minP, max: maxP };
            } catch (e) {
                console.warn('Failed to compute plane indices:', e);
                state.geom.planeId = null;
            }
        }

        // Fit the view once as soon as geometry is ready
        if (!state.viewFitted) {
            requestAnimationFrame(() => {
                if (!state.viewFitted && state.deckgl) {
                    window.pciSeq.rendering.autoFitView();
                    state.viewFitted = true;
                }
            });
        }
    });

    // Classes/probability streaming per iteration
    socket.on('classes_update_begin', (meta) => {
        console.log('=== BEGIN classes_update ===', meta);
        document.getElementById('loading-overlay').classList.add('hidden');
        state.stream = {
            iteration: meta.iteration,
            delta: meta.delta,
            total: meta.num_cells,
            received: 0,
            cell_classes: new Uint8Array(meta.num_cells),
            prob: new Float32Array(meta.num_cells)
        };
    });

    socket.on('classes_update_chunk', (chunk) => {
        if (!state.stream) return;
        const {start, end} = chunk;
        state.stream.cell_classes.set(chunk.cell_classes, start);
        state.stream.prob.set(chunk.prob, start);
        state.stream.received = end;
    });

    socket.on('classes_update_end', (msg) => {
        if (!state.stream) return;
        if (!state.geom.ready) {
            console.warn('Received classes without geometry; waiting.');
            return;
        }
        const N = state.stream.total;

        console.log(`=== CLASSES_UPDATE_END (iter ${state.stream.iteration}) ===`);
        console.log(`N=${N}, geom arrays length: centroids_x=${state.geom.centroids_x.length}, centroids_y=${state.geom.centroids_y.length}, radii=${state.geom.radii.length}`);

        const cells = new Array(N);
        for (let i = 0; i < N; i++) {
            cells[i] = {
                id: state.geom.cell_ids ? state.geom.cell_ids[i] : i,  // Use original_label from server, fallback to index
                x: state.geom.centroids_x[i],
                y: state.geom.centroids_y[i],
                radius: state.geom.radii[i],
                class: state.stream.cell_classes[i],
                prob: state.stream.prob[i]
            };
        }

        // Save OLD cells probability as previous BEFORE updating to new iteration
        if (state.cells && state.cells.length > 0) {
            state.previousProb = state.cells.map(cell => cell.prob);
            // Also save previous class assignments for highlight detection
            state.previousClass = state.cells.map(cell => cell.class);
        }

        state.cells = cells;
        state.numCells = N;
        state.iteration = state.stream.iteration;
        state.delta = state.stream.delta;
        state.stream = null;

        // Detect changed cells for changes view
        window.pciSeq.detection.detectChangedCells(false);

        // Detect cells that changed CLASS for highlighting
        window.pciSeq.detection.detectClassChanges();

        // Add to convergence chart history
        state.deltaHistory.push({ iteration: state.iteration, delta: state.delta });

        // Extend X-axis if needed
        if (state.iteration > state.chartXMax) {
            state.chartXMax = state.iteration + 2;
        }

        window.pciSeq.updateStatus();
        window.pciSeq.updateCellClassCounts();
        window.pciSeq.updateLegend();
        window.pciSeq.render();
        window.pciSeq.chart.renderConvergenceChart();
    });

    // Helper function for backward compatibility
    function applyClassesUpdate(iteration, delta, numCells, classes, prob) {
        if (!state.geom.ready) return;
        const N = numCells;
        const cells = new Array(N);
        for (let i = 0; i < N; i++) {
            cells[i] = {
                id: state.geom.cell_ids ? state.geom.cell_ids[i] : i,  // Use original_label from server, fallback to index
                x: state.geom.centroids_x[i],
                y: state.geom.centroids_y[i],
                radius: state.geom.radii[i],
                class: classes[i],
                prob: prob ? prob[i] : 1.0
            };
        }
        state.cells = cells;
        state.numCells = N;
        state.iteration = iteration;
        state.delta = delta;
        window.pciSeq.updateStatus();
        window.pciSeq.updateCellClassCounts();
        window.pciSeq.updateLegend();
        window.pciSeq.render();
    }

    // Diagnostics (check_cell) response listener is added in a separate commit

    // Export socket reference
    window.pciSeq.socket = socket;

})();
