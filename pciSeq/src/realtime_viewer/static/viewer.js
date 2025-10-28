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
    connected: false,
    deckgl: null,
    stream: null,  // holds buffers during chunked transfer
    geom: {
        ready: false,
        centroids_x: null,
        centroids_y: null,
        radii: null,
        stream: null
    }
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

    // Store class names mapping (index to name)
    if (meta.class_names) {
        state.cellClassNames = {};
        meta.class_names.forEach((name, idx) => {
            state.cellClassNames[idx] = name;
        });
        console.log('Loaded class names:', state.cellClassNames);
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
            x: state.geom.centroids_x[i],
            y: state.geom.centroids_y[i],
            radius: state.geom.radii[i],
            class: state.stream.cell_classes[i],
            confidence: state.stream.confidence[i]
        };
    }

    console.log(`Created ${cells.length} cells, sample cells (first 3):`);
    console.log(cells.slice(0, 3));

    state.cells = cells;
    state.numCells = N;
    state.iteration = state.stream.iteration;
    state.delta = state.stream.delta;
    state.stream = null;

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

    if (connected) {
        dot.classList.add('connected');
        status.textContent = 'Connected';
        status.style.color = '#00ff00';
    } else {
        dot.classList.remove('connected');
        status.textContent = 'Disconnected';
        status.style.color = '#ff0000';
    }
}

function updateStatus() {
    document.getElementById('iteration-value').textContent = state.iteration;
    document.getElementById('delta-value').textContent = state.delta.toFixed(6);
    document.getElementById('cells-value').textContent = state.numCells.toLocaleString();
}

function updateCellClassCounts() {
    if (!state.cells || state.cells.length === 0) return;

    // Count cells per class
    state.cellClassCounts = {};
    state.cells.forEach(cell => {
        const classIdx = cell.class;
        state.cellClassCounts[classIdx] = (state.cellClassCounts[classIdx] || 0) + 1;
    });
}

function updateLegend() {
    const legendItems = document.getElementById('legend-items');
    legendItems.innerHTML = '';

    // Sort classes by count (descending)
    const sortedClasses = Object.entries(state.cellClassCounts)
        .sort((a, b) => b[1] - a[1]);

    sortedClasses.forEach(([classIdx, count]) => {
        const item = document.createElement('div');
        item.className = 'legend-item';

        const colorBox = document.createElement('div');
        colorBox.className = 'legend-color';
        const color = state.cellClassColors[classIdx];
        colorBox.style.background = `rgb(${color[0]}, ${color[1]}, ${color[2]})`;

        const label = document.createElement('span');
        label.className = 'legend-label';
        // Use class name if available, otherwise fall back to index
        const className = state.cellClassNames[classIdx] || `Class ${classIdx}`;
        label.textContent = className;

        const countSpan = document.createElement('span');
        countSpan.className = 'legend-count';
        countSpan.textContent = count.toLocaleString();

        item.appendChild(colorBox);
        item.appendChild(label);
        item.appendChild(countSpan);
        legendItems.appendChild(item);
    });
}

function render() {
    if (!state.deckgl || state.cells.length === 0) {
        console.warn(`render() skipped: deckgl=${!!state.deckgl}, cells.length=${state.cells.length}`);
        return;
    }

    console.log(`=== RENDER (iter ${state.iteration}) ===`);
    console.log(`Rendering ${state.cells.length} cells`);
    console.log(`Sample cell data (first 3):`, state.cells.slice(0, 3));

    const {ScatterplotLayer} = deck;

    // Create scatterplot layer with cells
    const layer = new ScatterplotLayer({
        id: 'cells-layer',
        data: state.cells,
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
            getFillColor: [state.iteration]  // Update colors when iteration changes
        }
    });

    console.log(`Created layer with ${state.cells.length} data points`);

    // Update deck.gl with new layer
    state.deckgl.setProps({
        layers: [layer]
    });

    // Auto-fit view on first render
    if (state.iteration === 1) {
        console.log('Calling autoFitView because iteration === 1');
        autoFitView();
    }
}

function autoFitView() {
    if (state.cells.length === 0) return;

    console.log('=== AUTO FIT VIEW ===');

    // Calculate bounds
    let minX = Infinity, minY = Infinity;
    let maxX = -Infinity, maxY = -Infinity;

    state.cells.forEach(cell => {
        minX = Math.min(minX, cell.x);
        minY = Math.min(minY, cell.y);
        maxX = Math.max(maxX, cell.x);
        maxY = Math.max(maxY, cell.y);
    });

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

// Initialize when page loads
window.addEventListener('load', initialize);

// Handle window resize
window.addEventListener('resize', () => {
    if (state.deckgl) {
        state.deckgl.setProps({
            width: '100%',
            height: '100%'
        });
    }
});
