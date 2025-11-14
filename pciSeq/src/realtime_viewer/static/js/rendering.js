/**
 * Rendering Module
 * Handles deck.gl initialization, layer creation, and visualization
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Initialize deck.gl
    function initializeDeck() {
        const {DeckGL, ScatterplotLayer, OrthographicView} = deck;

        // Calculate initial view center from img_dim if available, otherwise use origin
        const defaultCenter = state.geom?.imgDim
            ? [state.geom.imgDim.w / 2, state.geom.imgDim.h / 2, 0]
            : [0, 0, 0];

        state.deckgl = new DeckGL({
            container: 'deck-container',
            views: [new OrthographicView({id: 'ortho', controller: true})],
            initialViewState: {
                target: defaultCenter,
                zoom: -1  // Start zoomed out
            },
            controller: true,
            layers: [],
            getTooltip: ({object}) => {
                if (object) {
                    const className = state.cellClassNames[object.class] || `Class ${object.class}`;
                    return {
                        html: `<div style="font-size: 12px;">
                            Cell Class: ${className}<br/>
                            Position: (${Math.round(object.x)}, ${Math.round(object.y)})<br/>
                            Prob: ${(object.prob * 100).toFixed(1)}%
                        </div>`,
                        style: {
                            backgroundColor: '#1b1b1b',
                            color: '#e5e5e5',
                            padding: '8px',
                            borderRadius: '4px'
                        }
                    };
                }
                return null;
            },
            onViewStateChange: ({viewState}) => {
                // Keep current view state so radius scaling can react to zoom
                state.viewState = viewState;
                return viewState;
            }
        });

        // Initialize stored view state for first render
        state.viewState = state.viewState || state.deckgl.props.initialViewState;

        render();
    }

    // Compute multiplicative scale so added on-screen radius ~ constant pixels
    function computeHighlightScale(baseRadius, opacity) {
        const zoom = (state.viewState && typeof state.viewState.zoom === 'number') ? state.viewState.zoom : -1;
        const screenScale = Math.pow(2, zoom);
        const pixelBoost = state.highlightPixelBoost || 12; // desired extra pixels at full opacity
        const denom = Math.max(baseRadius * screenScale, 1e-6);
        const ratio = Math.max(0.3, pixelBoost / denom); //clamp to 0.3

        // try {
        //     console.log(`highlight-scale debug: opacity=${Number(opacity).toFixed(3)}, ratio=${Number(ratio).toExponential(3)}`);
        // } catch (_) {
        //     // noop if formatting fails
        // }
        return 1.0 + (opacity * ratio);
    }

    // Main render function
    function render() {
        if (!state.deckgl || state.cells.length === 0) {
            console.warn(`render() skipped: deckgl=${!!state.deckgl}, cells.length=${state.cells.length}`);
            return;
        }

        // Filter cells based on visibility and view mode
        let visibleCells = state.cells.filter(cell => state.cellClassVisible[cell.class]);

        // Further filter by changes if in changes mode
        if (state.viewMode === 'changes') {
            visibleCells = visibleCells.filter(cell => state.changedCells.has(cell.id));
        }

        // Further filter by plane if enabled
        if (state.planeFilterEnabled && state.geom && state.geom.planeId && state.selectedPlane !== null) {
            visibleCells = visibleCells.filter(cell => state.geom.planeId[cell.id] === state.selectedPlane);
        }

        console.log(`=== RENDER (iter ${state.iteration}) ===`);
        const modeText = state.viewMode === 'changes' ? 'changes mode' : 'all cells mode';
        console.log(`Rendering ${visibleCells.length}/${state.cells.length} cells (${modeText})`);
        console.log(`Sample cell data (first 3):`, visibleCells.slice(0, 3));

        const {ScatterplotLayer} = deck;

        // Get current timestamp for fade calculations
        const now = Date.now();

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
            getRadius: d => {
                const baseRadius = state.geom && state.geom.mcr ? state.geom.mcr : d.radius;

                // Check if this cell changed class and should be scaled up
                const changeTime = state.classChangedCells.get(d.id);
                if (changeTime) {
                    const elapsed = now - changeTime;
                    const fadeDuration = state.highlightFadeDuration;
                    const fadeProgress = elapsed / fadeDuration; // 0 to 1
                    const opacity = Math.max(0, 1 - fadeProgress); // 1 to 0

                    if (opacity > 0) {
                        const scale = computeHighlightScale(baseRadius, opacity);
                        return baseRadius * scale;
                    }
                }

                return baseRadius;
            },
            getFillColor: d => {
                const color = state.cellClassColors[d.class] || [128, 128, 128];
                // Make all cells clearly visible: clamp alpha to [0.7, 1.0]
                const alpha = Math.round((0.7 + d.prob * 0.3) * 255);
                return [color[0], color[1], color[2], alpha];
            },
            getLineColor: [255, 255, 255, 60],
            updateTriggers: {
                getFillColor: [state.iteration],
                getRadius: [state.iteration, state.classChangedCells.size, state.viewState ? state.viewState.zoom : 0],
                data: [Object.values(state.cellClassVisible)]
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

    // Auto-fit view to show all cells
    function autoFitView() {
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

    // Update legend with current cell class counts
    window.pciSeq.updateLegend = function() {
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

            // Subtle eye/eye-off icon (visual cue)
            const eyeWrap = document.createElement('span');
            eyeWrap.className = 'chip-eye';
            const eyeOpenSvg = `
                <svg viewBox="0 0 24 24" aria-hidden="true">
                  <path d="M12 5c-7 0-11 7-11 7s4 7 11 7 11-7 11-7-4-7-11-7z" fill="none" stroke="currentColor" stroke-width="1.5"/>
                  <circle cx="12" cy="12" r="3" fill="currentColor"/>
                </svg>`;
            const eyeOffSvg = `
                <svg viewBox="0 0 24 24" aria-hidden="true">
                  <path d="M12 5c-7 0-11 7-11 7s4 7 11 7 11-7 11-7-4-7-11-7z" fill="none" stroke="currentColor" stroke-width="1.5"/>
                  <circle cx="12" cy="12" r="3" fill="currentColor"/>
                  <line x1="4" y1="4" x2="20" y2="20" stroke="currentColor" stroke-width="1.5"/>
                </svg>`;
            eyeWrap.innerHTML = isVisible ? eyeOpenSvg : eyeOffSvg;

            // Assemble
            item.appendChild(colorBox);
            item.appendChild(label);
            item.appendChild(countSpan);
            item.appendChild(eyeWrap);

            // Click toggles visibility
            item.addEventListener('click', () => toggleClassVisibility(classIdx));

            legendItems.appendChild(item);
        });
    };

    // Toggle class visibility
    function toggleClassVisibility(classIdx) {
        state.cellClassVisible[classIdx] = !state.cellClassVisible[classIdx];
        window.pciSeq.updateLegend();
        window.pciSeq.updateStatus();
        render();
    }

    // Show all classes
    function showAllClasses() {
        Object.keys(state.cellClassVisible).forEach(classIdx => {
            state.cellClassVisible[classIdx] = true;
        });
        window.pciSeq.updateLegend();
        window.pciSeq.updateStatus();
        render();
    }

    // Hide all classes
    function hideAllClasses() {
        Object.keys(state.cellClassVisible).forEach(classIdx => {
            state.cellClassVisible[classIdx] = false;
        });
        window.pciSeq.updateLegend();
        window.pciSeq.updateStatus();
        render();
    }

    // Export functions
    window.pciSeq.render = render;
    window.pciSeq.rendering = {
        initializeDeck: initializeDeck,
        autoFitView: autoFitView,
        showAllClasses: showAllClasses,
        hideAllClasses: hideAllClasses
    };

})();
