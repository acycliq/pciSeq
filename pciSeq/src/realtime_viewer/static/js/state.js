/**
 * State Management Module
 * Centralized state for the pciSeq Real-Time Viewer
 */

(function() {
    'use strict';

    // Initialize global namespace
    window.pciSeq = window.pciSeq || {};

    // Application state
    const state = {
        cells: [],  // Array of {x, y, radius, class, prob}
        numCells: 0,
        iteration: 0,
        delta: 0,
        cellClassColors: {},
        cellClassCounts: {},
        cellClassNames: {},  // Maps class index to class name
        cellClassVisible: {},  // Maps class index to visibility (true/false)
        pendingColorScheme: null,  // Store color scheme loaded before class names arrive

        // Changes view state
        viewMode: 'all',  // 'all' or 'changes'
        changeThreshold: 1.0,  // Percentage threshold for significant change
        previousProb: null,  // Store previous iteration's probability values
        changedCells: new Set(),  // Set of cell indices that changed

        // Class change highlighting (fade-out effect)
        previousClass: null,  // Store previous iteration's class assignments
        classChangedCells: new Map(),  // Map of cell_id -> timestamp for cells that changed class
        animationFrameId: null,  // Track animation loop
        highlightFadeDuration: 2000,  // Fade duration in milliseconds
        highlightColor: [0, 217, 255],  // Soft cyan RGB

        // Convergence chart
        cellCallTolerance: 0.2,  // Convergence threshold (from server)
        deltaHistory: [],  // Array of {iteration, delta} objects
        chartXMax: 10,  // Current X-axis max (dynamically extends)

        connected: false,
        deckgl: null,
        stream: null,  // holds buffers during chunked transfer
        geom: {
            ready: false,
            centroids_x: null,
            centroids_y: null,
            centroids_z: null,
            radii: null,
            stream: null,
            mcr: null,
            is3D: false,
            voxelSize: null,
            planeId: null
        },
        // Has the view been auto-fitted once?
        viewFitted: false,
        // Plane filter (3D)
        planeFilterEnabled: false,
        selectedPlane: null,
        planeRange: null,
        // Quick legend filter text (appears on demand)
        legendFilter: ''
    };

    // Export state and helper functions
    window.pciSeq.state = state;

    // State update functions
    window.pciSeq.updateConnectionStatus = function(connected) {
        state.connected = connected;

        const dot = document.getElementById('connection-dot');
        const status = document.getElementById('connection-status');
        const hudDot = document.getElementById('hud-connection-dot');

        // Read CSS variables for consistent colors
        const cssVars = getComputedStyle(document.documentElement);
        const accent = (cssVars.getPropertyValue('--accent') || '#22c55e').trim();
        const danger = (cssVars.getPropertyValue('--danger') || '#ef4444').trim();

        if (connected) {
            dot.classList.add('connected');
            status.textContent = 'Connected';
            status.style.color = accent;
            if (hudDot) {
                hudDot.classList.add('connected');
                hudDot.style.background = accent;
            }
        } else {
            dot.classList.remove('connected');
            status.textContent = 'Disconnected';
            status.style.color = danger;
            if (hudDot) {
                hudDot.classList.remove('connected');
                hudDot.style.background = danger;
            }
        }
    };

    window.pciSeq.updateStatus = function() {
        document.getElementById('iteration-value').textContent = state.iteration;
        document.getElementById('delta-value').textContent = state.delta.toFixed(3);
        document.getElementById('cells-value').textContent = state.numCells.toLocaleString();

        // Update the minimal HUD (only cells count now)
        const hCells = document.getElementById('hud-cells-value');
        if (hCells) hCells.textContent = state.numCells.toLocaleString();
    };

    window.pciSeq.updateCellClassCounts = function() {
        if (!state.cells || state.cells.length === 0) return;

        // Count cells per class
        state.cellClassCounts = {};
        state.cells.forEach(cell => {
            const classIdx = cell.class;
            state.cellClassCounts[classIdx] = (state.cellClassCounts[classIdx] || 0) + 1;

            // Initialize visibility to true for new classes
            if (!(classIdx in state.cellClassVisible)) {
                state.cellClassVisible[classIdx] = true;
            }
        });
    };

    window.pciSeq.updateChangesCount = function() {
        const countEl = document.getElementById('changes-count');
        if (countEl) {
            countEl.textContent = state.changedCells.size.toLocaleString();
        }
    };

    // Simple user-facing notice (non-technical) shown as a small banner
    window.pciSeq.showUserNotice = function(message) {
        let el = document.getElementById('user-notice');
        if (!el) {
            el = document.createElement('div');
            el.id = 'user-notice';
            el.style.position = 'fixed';
            el.style.bottom = '20px';
            el.style.left = '20px';
            el.style.maxWidth = '420px';
            el.style.padding = '10px 12px';
            el.style.background = '#2a2a2a';
            el.style.color = '#e0e0e0';
            el.style.border = '1px solid #00ff00';
            el.style.borderRadius = '4px';
            el.style.boxShadow = '0 2px 8px rgba(0,0,0,0.4)';
            el.style.fontSize = '12px';
            el.style.zIndex = 2000;
            document.body.appendChild(el);
        }
        el.textContent = message;
        el.style.display = 'block';
        clearTimeout(el._hideTimer);
        el._hideTimer = setTimeout(() => { el.style.display = 'none'; }, 7000);
    };

})();