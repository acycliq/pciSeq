/**
 * pciSeq Real-Time Algorithm Viewer - Spatial View with deck.gl
 * Visualizes cells at their actual spatial positions as they update during VarBayes execution
 */

// State
const state = {
    cells: [],  // Array of {x, y, radius, class, confidence}
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
    previousConfidence: null,  // Store previous iteration's confidence values
    changedCells: new Set(),  // Set of cell indices that changed

    connected: false,
    deckgl: null,
    stream: null,  // holds buffers during chunked transfer
    geom: {
        ready: false,
        centroids_x: null,
        centroids_y: null,
        radii: null,
        stream: null
    },
    // Has the view been auto-fitted once?
    viewFitted: false,
    // Quick legend filter text (appears on demand)
    legendFilter: ''
};

// Socket.IO connection
const socket = io();

// Event handlers
socket.on('connect', () => {
    console.log('Connected to pciSeq server');
    state.connected = true;
    updateConnectionStatus(true);
});

socket.on('disconnect', (reason) => {
    console.log('Disconnected:', reason);
    state.connected = false;
    updateConnectionStatus(false);
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
        applyClassesUpdate(data.iteration, data.delta, data.num_cells, data.cell_classes, data.confidence);
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
        centroids_x: new Float32Array(meta.num_cells),
        centroids_y: new Float32Array(meta.num_cells),
        radii: new Float32Array(meta.num_cells)
    };

    // Reset change tracking at the start of a new run/session.
    // This prevents comparing the new run against old data after reconnects or restarts.
    state.previousConfidence = null;   // No previous step yet for the new run
    state.changedCells.clear();        // Start with an empty set of changed cells
    updateChangesCount();              // Reflect the reset in the UI counter
    // Hide any previous informational banner, if present
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
            applyColorScheme(state.pendingColorScheme);
            state.pendingColorScheme = null;
        }
    }
});

socket.on('geometry_init_chunk', (chunk) => {
    if (!state.geom.stream) return;
    const {start, end} = chunk;
    state.geom.stream.centroids_x.set(chunk.centroids_x, start);
    state.geom.stream.centroids_y.set(chunk.centroids_y, start);
    state.geom.stream.radii.set(chunk.radii, start);
    state.geom.stream.received = end;
});

socket.on('geometry_init_end', () => {
    console.log('=== END geometry_init ===');
    const gs = state.geom.stream;
    state.geom.centroids_x = gs.centroids_x;
    state.geom.centroids_y = gs.centroids_y;
    state.geom.radii = gs.radii;
    state.geom.ready = true;
    state.geom.stream = null;

    // Fit the view once as soon as geometry is ready, so the
    // initial chart is centered and justified before the first update.
    if (!state.viewFitted) {
        // Defer to the next frame to ensure deck.gl is initialized
        requestAnimationFrame(() => {
            if (!state.viewFitted && state.deckgl) {
                autoFitView();
                state.viewFitted = true;
            }
        });
    }
});

// Classes/confidence streaming per iteration
socket.on('classes_update_begin', (meta) => {
    console.log('=== BEGIN classes_update ===', meta);
    document.getElementById('loading-overlay').classList.add('hidden');
    state.stream = {
        iteration: meta.iteration,
        delta: meta.delta,
        total: meta.num_cells,
        received: 0,
        cell_classes: new Uint8Array(meta.num_cells),
        confidence: new Float32Array(meta.num_cells)
    };
});

socket.on('classes_update_chunk', (chunk) => {
    if (!state.stream) return;
    const {start, end} = chunk;
    state.stream.cell_classes.set(chunk.cell_classes, start);
    state.stream.confidence.set(chunk.confidence, start);
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
    console.log(`Sample centroids (first 3): x=[${state.geom.centroids_x[0]}, ${state.geom.centroids_x[1]}, ${state.geom.centroids_x[2]}], y=[${state.geom.centroids_y[0]}, ${state.geom.centroids_y[1]}, ${state.geom.centroids_y[2]}]`);
    console.log(`Sample radii (first 3): [${state.geom.radii[0]}, ${state.geom.radii[1]}, ${state.geom.radii[2]}]`);
    console.log(`Sample classes (first 3): [${state.stream.cell_classes[0]}, ${state.stream.cell_classes[1]}, ${state.stream.cell_classes[2]}]`);

    const cells = new Array(N);
    for (let i = 0; i < N; i++) {
        cells[i] = {
            id: i,  // Stable cell ID for tracking across filtering operations
            x: state.geom.centroids_x[i],
            y: state.geom.centroids_y[i],
            radius: state.geom.radii[i],
            class: state.stream.cell_classes[i],
            confidence: state.stream.confidence[i]
        };
    }

    console.log(`Created ${cells.length} cells, sample cells (first 3):`);
    console.log(cells.slice(0, 3));

    // Save OLD cells confidence as previous BEFORE updating to new iteration
    if (state.cells && state.cells.length > 0) {
        state.previousConfidence = state.cells.map(cell => cell.confidence);
    }

    state.cells = cells;
    state.numCells = N;
    state.iteration = state.stream.iteration;
    state.delta = state.stream.delta;
    state.stream = null;

    // Detect changed cells for changes view (pass false since we already saved previous)
    detectChangedCells(false);

    updateStatus();
    updateCellClassCounts();
    updateLegend();
    render();
});

function applyClassesUpdate(iteration, delta, numCells, classes, confidence) {
    if (!state.geom.ready) return;
    const N = numCells;
    const cells = new Array(N);
    for (let i = 0; i < N; i++) {
        cells[i] = {
            id: i,  // Stable cell ID for tracking across filtering operations
            x: state.geom.centroids_x[i],
            y: state.geom.centroids_y[i],
            radius: state.geom.radii[i],
            class: classes[i],
            confidence: confidence ? confidence[i] : 1.0
        };
    }
    state.cells = cells;
    state.numCells = N;
    state.iteration = iteration;
    state.delta = delta;
    updateStatus();
    updateCellClassCounts();
    updateLegend();
    render();
}

// Initialize visualization
function initialize() {
    // Generate color palette for cell classes
    generateColorPalette();

    // Initialize deck.gl
    initializeDeck();
}

function initializeDeck() {
    const {DeckGL, ScatterplotLayer, OrthographicView} = deck;

    state.deckgl = new DeckGL({
        container: 'deck-container',
        views: [new OrthographicView({id: 'ortho', controller: true})],
        initialViewState: {
            target: [3200, 2200, 0],  // Center of typical image
            zoom: -1  // Start zoomed out
        },
        controller: true,
        layers: [],
        getTooltip: ({object}) => {
            if (object) {
                const className = state.cellClassNames[object.class] || `Class ${object.class}`;
                return {
                    html: `<div style="font-family: monospace; font-size: 12px;">
                        Cell Class: ${className}<br/>
                        Position: (${Math.round(object.x)}, ${Math.round(object.y)})<br/>
                        Confidence: ${(object.confidence * 100).toFixed(1)}%
                    </div>`,
                    style: {
                        backgroundColor: '#2a2a2a',
                        color: '#00ff00',
                        padding: '8px',
                        borderRadius: '4px'
                    }
                };
            }
            return null;
        },
        onViewStateChange: ({viewState}) => {
            // Optional: store view state if needed
            return viewState;
        }
    });

    render();
}

function generateColorPalette() {
    // Generate distinct colors for up to 65 classes
    const numColors = 65;
    const hueStep = 360 / numColors;

    for (let i = 0; i < numColors; i++) {
        const hue = (i * hueStep) % 360;
        const saturation = 70 + (i % 3) * 10; // Vary saturation slightly
        const lightness = 50 + (i % 2) * 10;  // Vary lightness slightly
        state.cellClassColors[i] = hslToRgb(hue, saturation, lightness);
    }
}

function hslToRgb(h, s, l) {
    s /= 100;
    l /= 100;
    const k = n => (n + h / 30) % 12;
    const a = s * Math.min(l, 1 - l);
    const f = n => l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
    return [
        Math.round(255 * f(0)),
        Math.round(255 * f(8)),
        Math.round(255 * f(4))
    ];
}

function updateConnectionStatus(connected) {
    const dot = document.getElementById('connection-dot');
    const status = document.getElementById('connection-status');
    const hudDot = document.getElementById('hud-connection-dot');

    if (connected) {
        dot.classList.add('connected');
        status.textContent = 'Connected';
        status.style.color = '#00ff00';
        if (hudDot) hudDot.classList.add('connected');
    } else {
        dot.classList.remove('connected');
        status.textContent = 'Disconnected';
        status.style.color = '#ff0000';
        if (hudDot) hudDot.classList.remove('connected');
    }
}

function updateStatus() {
    document.getElementById('iteration-value').textContent = state.iteration;
    document.getElementById('delta-value').textContent = state.delta.toFixed(6);
    document.getElementById('cells-value').textContent = state.numCells.toLocaleString();

    // Also update the minimal HUD if present
    const hIter = document.getElementById('hud-iteration-value');
    const hDelta = document.getElementById('hud-delta-value');
    const hCells = document.getElementById('hud-cells-value');
    if (hIter) hIter.textContent = state.iteration;
    if (hDelta) hDelta.textContent = state.delta.toFixed(6);
    if (hCells) hCells.textContent = state.numCells.toLocaleString();
}

function updateCellClassCounts() {
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
}

function updateLegend() {
    const legendItems = document.getElementById('legend-items');
    legendItems.innerHTML = '';

    // Filter and sort classes
    let entries = Object.entries(state.cellClassCounts);
    const filterText = (state.legendFilter || '').trim().toLowerCase();
    if (filterText) {
        entries = entries.filter(([idx]) => {
            const name = state.cellClassNames[idx] || `Class ${idx}`;
            return String(name).toLowerCase().includes(filterText);
        });
    }
    const sortedClasses = entries.sort((a, b) => b[1] - a[1]);

    sortedClasses.forEach(([classIdx, count]) => {
        const item = document.createElement('div');
        item.className = 'chip';
        const isVisible = state.cellClassVisible[classIdx];
        if (!isVisible) item.classList.add('dim');

        // Color swatch
        const colorBox = document.createElement('div');
        colorBox.className = 'legend-color';
        const color = state.cellClassColors[classIdx] || [128,128,128];
        colorBox.style.background = `rgb(${color[0]}, ${color[1]}, ${color[2]})`;

        // Label
        const label = document.createElement('span');
        label.className = 'legend-label';
        const className = state.cellClassNames[classIdx] || `Class ${classIdx}`;
        label.textContent = className;

        // Count (muted, with tooltip)
        const countSpan = document.createElement('span');
        countSpan.className = 'legend-count';
        countSpan.textContent = count.toLocaleString();
        item.title = `${className}: ${count.toLocaleString()} cells`;

        // Assemble
        item.appendChild(colorBox);
        item.appendChild(label);
        item.appendChild(countSpan);

        // Click toggles visibility
        item.addEventListener('click', () => toggleClassVisibility(classIdx));

        legendItems.appendChild(item);
    });
}

function toggleClassVisibility(classIdx) {
    state.cellClassVisible[classIdx] = !state.cellClassVisible[classIdx];
    updateLegend();
    updateStatus();
    render();
}

function showAllClasses() {
    Object.keys(state.cellClassVisible).forEach(classIdx => {
        state.cellClassVisible[classIdx] = true;
    });
    updateLegend();
    updateStatus();
    render();
}

function hideAllClasses() {
    Object.keys(state.cellClassVisible).forEach(classIdx => {
        state.cellClassVisible[classIdx] = false;
    });
    updateLegend();
    updateStatus();
    render();
}

function detectChangedCells(updatePrevious = true) {
    // Clear changed cells set
    state.changedCells.clear();

    // Safety guard: if the number of cells changed between steps,
    // pause change detection for this step and reset the baseline.
    if (state.previousConfidence && state.previousConfidence.length !== state.cells.length) {
        const prevN = state.previousConfidence.length;
        const currN = state.cells.length;
        const msg = `Paused the "Changes" view for this step because the number of cells changed (${prevN} → ${currN}). This can happen after reconnecting or restarting. The view will resume automatically on the next step.`;
        console.warn(msg);
        showUserNotice(msg);
        // Reset baseline to current values and skip diffing for this step
        state.previousConfidence = state.cells.map(cell => cell.confidence);
        updateChangesCount();
        return;
    }

    // On first iteration (no previous data to compare), all cells are "changed"
    if (!state.previousConfidence) {
        for (let i = 0; i < state.cells.length; i++) {
            state.changedCells.add(i);
        }
        console.log(`First iteration (no previous data): marked all ${state.cells.length} cells as changed`);
    } else {
        // Compare with previous iteration
        const threshold = state.changeThreshold / 100; // Convert percentage to decimal

        for (let i = 0; i < state.cells.length; i++) {
            const currentConf = state.cells[i].confidence;
            const prevConf = state.previousConfidence[i];

            // Calculate absolute change in confidence
            const change = Math.abs(currentConf - prevConf);

            // Mark as changed if exceeds threshold
            if (change >= threshold) {
                state.changedCells.add(i);
            }
        }

        console.log(`Detected ${state.changedCells.size}/${state.cells.length} changed cells (threshold: ${state.changeThreshold}%)`);
    }

    // Only store current confidence when new iteration arrives, not when threshold changes
    if (updatePrevious) {
        state.previousConfidence = state.cells.map(cell => cell.confidence);
    }

    // Update changes count display
    updateChangesCount();
}

// Simple user-facing notice (non-technical) shown as a small banner
function showUserNotice(message) {
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
}

function updateChangesCount() {
    const countEl = document.getElementById('changes-count');
    if (countEl) {
        countEl.textContent = state.changedCells.size.toLocaleString();
    }
}

function render() {
    if (!state.deckgl || state.cells.length === 0) {
        console.warn(`render() skipped: deckgl=${!!state.deckgl}, cells.length=${state.cells.length}`);
        return;
    }

    // Filter cells based on visibility and view mode
    let visibleCells = state.cells.filter(cell => state.cellClassVisible[cell.class]);

    // Further filter by changes if in changes mode
    if (state.viewMode === 'changes') {
        // Use original cell id for membership test, not filtered index
        visibleCells = visibleCells.filter(cell => state.changedCells.has(cell.id));
    }

    console.log(`=== RENDER (iter ${state.iteration}) ===`);
    const modeText = state.viewMode === 'changes' ? 'changes mode' : 'all cells mode';
    console.log(`Rendering ${visibleCells.length}/${state.cells.length} cells (${modeText})`);
    console.log(`Sample cell data (first 3):`, visibleCells.slice(0, 3));

    const {ScatterplotLayer} = deck;

    // Create scatterplot layer with visible cells only
    const layer = new ScatterplotLayer({
        id: 'cells-layer',
        data: visibleCells,
        pickable: true,
        opacity: 1.0,
        stroked: true,
        filled: true,
        radiusScale: 1,
        radiusMinPixels: 2,
        radiusMaxPixels: 100,
        lineWidthMinPixels: 1,
        getPosition: d => [d.x, d.y],
        getRadius: d => d.radius,
        getFillColor: d => {
            const color = state.cellClassColors[d.class] || [128, 128, 128];
            // Make all cells clearly visible: clamp alpha to [0.7, 1.0]
            const alpha = Math.round((0.7 + d.confidence * 0.3) * 255);
            return [color[0], color[1], color[2], alpha];
        },
        getLineColor: [255, 255, 255, 60],
        updateTriggers: {
            getFillColor: [state.iteration],  // Update colors when iteration changes
            data: [Object.values(state.cellClassVisible)]  // Update when visibility changes
        }
    });

    console.log(`Created layer with ${visibleCells.length} visible data points`);

    // Update deck.gl with new layer
    state.deckgl.setProps({
        layers: [layer]
    });

    // Auto-fit view on the very first render if it hasn't been fitted yet
    if (!state.viewFitted) {
        console.log('Auto-fitting view on first render');
        autoFitView();
        state.viewFitted = true;
    }
}

function autoFitView() {
    // Compute bounds either from current cells or, if not set yet,
    // from the geometry centroids received during initialization.
    const useCells = state.cells && state.cells.length > 0;
    const useGeom = !useCells && state.geom && state.geom.ready && state.geom.centroids_x && state.geom.centroids_x.length > 0;
    if (!useCells && !useGeom) return;

    console.log('=== AUTO FIT VIEW ===');

    // Calculate bounds
    let minX = Infinity, minY = Infinity;
    let maxX = -Infinity, maxY = -Infinity;

    if (useCells) {
        state.cells.forEach(cell => {
            minX = Math.min(minX, cell.x);
            minY = Math.min(minY, cell.y);
            maxX = Math.max(maxX, cell.x);
            maxY = Math.max(maxY, cell.y);
        });
    } else if (useGeom) {
        const xs = state.geom.centroids_x;
        const ys = state.geom.centroids_y;
        for (let i = 0; i < xs.length; i++) {
            const x = xs[i];
            const y = ys[i];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
        }
    }

    console.log(`Bounds: minX=${minX}, maxX=${maxX}, minY=${minY}, maxY=${maxY}`);

    // Add 10% padding
    const width = maxX - minX;
    const height = maxY - minY;
    const paddedMinX = minX - width * 0.1;
    const paddedMinY = minY - height * 0.1;
    const paddedMaxX = maxX + width * 0.1;
    const paddedMaxY = maxY + height * 0.1;

    // Calculate center
    const centerX = (paddedMinX + paddedMaxX) / 2;
    const centerY = (paddedMinY + paddedMaxY) / 2;

    // Calculate zoom to fit
    const container = document.getElementById('deck-container');
    const containerWidth = container.clientWidth;
    const containerHeight = container.clientHeight;

    const dataWidth = paddedMaxX - paddedMinX;
    const dataHeight = paddedMaxY - paddedMinY;

    const zoomX = Math.log2(containerWidth / dataWidth);
    const zoomY = Math.log2(containerHeight / dataHeight);
    const zoom = Math.min(zoomX, zoomY);

    console.log(`Calculated view: center=[${centerX}, ${centerY}], zoom=${zoom}`);

    // Update view to fit all cells
    state.deckgl.setProps({
        initialViewState: {
            target: [centerX, centerY, 0],
            zoom: zoom,
            transitionDuration: 1000
        }
    });
}

// Custom color scheme import
function hexToRgb(hex) {
    // Remove # if present
    hex = hex.replace(/^#/, '');

    // Parse hex values
    const bigint = parseInt(hex, 16);
    const r = (bigint >> 16) & 255;
    const g = (bigint >> 8) & 255;
    const b = bigint & 255;

    return [r, g, b];
}

function applyColorScheme(colorScheme) {
    let appliedCount = 0;
    const notFoundClasses = [];

    // Create reverse mapping: class name -> class index
    const nameToIndex = {};
    Object.entries(state.cellClassNames).forEach(([idx, name]) => {
        nameToIndex[name] = parseInt(idx);
    });

    // Apply custom colors
    Object.entries(colorScheme).forEach(([className, hexColor]) => {
        const classIdx = nameToIndex[className];

        if (classIdx !== undefined) {
            try {
                const rgb = hexToRgb(hexColor);
                state.cellClassColors[classIdx] = rgb;
                appliedCount++;
            } catch (err) {
                console.warn(`Invalid color format for ${className}: ${hexColor}`);
            }
        } else {
            notFoundClasses.push(className);
        }
    });

    if (notFoundClasses.length > 0) {
        console.warn(`Classes not found in data: ${notFoundClasses.join(', ')}`);
    }

    return { appliedCount, notFoundClasses };
}

function loadCustomColors(colorScheme) {
    const statusEl = document.getElementById('file-status');

    // Check if class names have been loaded yet
    const hasClassNames = Object.keys(state.cellClassNames).length > 0;

    if (!hasClassNames) {
        // Store color scheme for later application
        state.pendingColorScheme = colorScheme;
        statusEl.textContent = `Color scheme loaded (${Object.keys(colorScheme).length} classes). Will apply when algorithm starts.`;
        statusEl.className = 'file-status success';

        // Clear status after 5 seconds
        setTimeout(() => {
            statusEl.textContent = '';
            statusEl.className = 'file-status';
        }, 5000);
        return;
    }

    // Apply colors immediately
    const { appliedCount, notFoundClasses } = applyColorScheme(colorScheme);

    // Update UI
    if (appliedCount > 0) {
        statusEl.textContent = `Applied ${appliedCount} custom colors`;
        statusEl.className = 'file-status success';

        // Refresh legend and visualization
        updateLegend();
        render();
    } else {
        statusEl.textContent = 'No matching classes found';
        statusEl.className = 'file-status error';
    }

    // Clear status after 5 seconds
    setTimeout(() => {
        statusEl.textContent = '';
        statusEl.className = 'file-status';
    }, 5000);
}

function handleColorFileUpload(event) {
    const file = event.target.files[0];
    const statusEl = document.getElementById('file-status');

    if (!file) return;

    statusEl.textContent = 'Loading...';
    statusEl.className = 'file-status';

    const reader = new FileReader();

    reader.onload = (e) => {
        try {
            const colorScheme = JSON.parse(e.target.result);

            // Validate JSON structure
            if (typeof colorScheme !== 'object' || Array.isArray(colorScheme)) {
                throw new Error('Invalid JSON format. Expected object with class_name: hex_color pairs');
            }

            loadCustomColors(colorScheme);
        } catch (err) {
            statusEl.textContent = `Error: ${err.message}`;
            statusEl.className = 'file-status error';
            console.error('Failed to load color scheme:', err);
        }
    };

    reader.onerror = () => {
        statusEl.textContent = 'Failed to read file';
        statusEl.className = 'file-status error';
    };

    reader.readAsText(file);

    // Reset file input so same file can be uploaded again
    event.target.value = '';
}

// Initialize when page loads
window.addEventListener('load', () => {
    initialize();

    // Setup color file upload handler
    const fileInput = document.getElementById('color-file-input');
    if (fileInput) {
        fileInput.addEventListener('change', handleColorFileUpload);
    }

    // Setup Show/Hide All buttons
    const showAllBtn = document.getElementById('show-all-btn');
    const hideAllBtn = document.getElementById('hide-all-btn');

    if (showAllBtn) {
        showAllBtn.addEventListener('click', showAllClasses);
    }

    if (hideAllBtn) {
        hideAllBtn.addEventListener('click', hideAllClasses);
    }

    // Setup compact Changes View segmented control
    const modeAllBtn = document.getElementById('mode-all');
    const modeChangesBtn = document.getElementById('mode-changes');
    const thresholdSection = document.getElementById('threshold-section');
    const changesInfo = document.getElementById('changes-info');

    function setViewMode(mode) {
        if (!mode || (mode !== 'all' && mode !== 'changes')) return;
        state.viewMode = mode;
        if (modeAllBtn && modeChangesBtn) {
            if (mode === 'all') {
                modeAllBtn.classList.add('active');
                modeChangesBtn.classList.remove('active');
            } else {
                modeAllBtn.classList.remove('active');
                modeChangesBtn.classList.add('active');
            }
        }
        // Show/hide threshold and changes count based on mode
        if (thresholdSection) thresholdSection.classList.toggle('hidden', mode !== 'changes');
        if (changesInfo) changesInfo.classList.toggle('hidden', mode !== 'changes');

        // If switching to changes, ensure changed set is computed with current threshold
        if (mode === 'changes' && state.cells.length > 0) {
            detectChangedCells(false);
        }
        render();
    }

    if (modeAllBtn) modeAllBtn.addEventListener('click', () => setViewMode('all'));
    if (modeChangesBtn) modeChangesBtn.addEventListener('click', () => setViewMode('changes'));
    // Initialize UI to current state.viewMode
    setViewMode(state.viewMode);

    const thresholdSlider = document.getElementById('threshold-slider');
    const thresholdValueDisplay = document.getElementById('threshold-value');

    if (thresholdSlider && thresholdValueDisplay) {
        thresholdSlider.addEventListener('input', (e) => {
            state.changeThreshold = parseFloat(e.target.value);
            thresholdValueDisplay.textContent = state.changeThreshold.toFixed(1);

            // Recalculate changed cells with new threshold, but DON'T update previous confidence
            if (state.cells.length > 0 && state.previousConfidence) {
                detectChangedCells(false);  // false = don't update previousConfidence
                render();
            }
        });
    }

    // Quick legend filter: toggle with '/', filter on input, hide on Esc/blur
    const legendFilterContainer = document.getElementById('legend-filter-container');
    const legendFilterInput = document.getElementById('legend-filter-input');

    // Shortcut: '/' focuses the filter input
    window.addEventListener('keydown', (e) => {
        if (e.key === '/' && !e.ctrlKey && !e.metaKey && !e.altKey) {
            const tag = document.activeElement && document.activeElement.tagName;
            if (tag === 'INPUT' || tag === 'TEXTAREA') return;
            e.preventDefault();
            if (legendFilterInput) {
                legendFilterInput.focus();
                legendFilterInput.select();
            }
        }
    });

    if (legendFilterInput) {
        legendFilterInput.addEventListener('input', () => {
            state.legendFilter = legendFilterInput.value;
            updateLegend();
        });
        legendFilterInput.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                legendFilterInput.value = '';
                state.legendFilter = '';
                updateLegend();
            }
        });
    }
});

// Handle window resize
window.addEventListener('resize', () => {
    if (state.deckgl) {
        state.deckgl.setProps({
            width: '100%',
            height: '100%'
        });
    }
});
