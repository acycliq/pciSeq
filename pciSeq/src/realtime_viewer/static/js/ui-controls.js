/**
 * UI Controls Module
 * Handles user interface interactions, buttons, sliders, and keyboard shortcuts
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Initialize UI controls
    function initializeControls() {
        setupColorFileUpload();
        setupShowHideButtons();
        setupViewModeControls();
        setupPlaneControls();
        setupThresholdSlider();
        setupLegendFilter();
        setupWindowResize();
    }

    // Setup color file upload handler
    function setupColorFileUpload() {
        const fileInput = document.getElementById('color-file-input');
        if (fileInput) {
            fileInput.addEventListener('change', window.pciSeq.colors.handleColorFileUpload);
        }
    }

    // Setup Show/Hide All buttons
    function setupShowHideButtons() {
        const showAllBtn = document.getElementById('show-all-btn');
        const hideAllBtn = document.getElementById('hide-all-btn');

        if (showAllBtn) {
            showAllBtn.addEventListener('click', window.pciSeq.rendering.showAllClasses);
        }

        if (hideAllBtn) {
            hideAllBtn.addEventListener('click', window.pciSeq.rendering.hideAllClasses);
        }
    }

    // Setup view mode controls (All / Changes)
    function setupViewModeControls() {
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
                window.pciSeq.detection.detectChangedCells(false);
            }
            window.pciSeq.render();
        }

        if (modeAllBtn) modeAllBtn.addEventListener('click', () => setViewMode('all'));
        if (modeChangesBtn) modeChangesBtn.addEventListener('click', () => setViewMode('changes'));
        // Initialize UI to current state.viewMode
        setViewMode(state.viewMode);
    }

    // Setup plane controls for 3D filtering
    function setupPlaneControls() {
        const planeToggleBtn = document.getElementById('plane-toggle');
        const planeControls = document.getElementById('plane-controls');
        const planeSlider = document.getElementById('planeSlider');
        const planeLabel = document.getElementById('planeLabel');

        function updatePlaneFooterVisibility() {
            if (!planeControls) return;
            const show = state.geom.is3D && state.planeFilterEnabled;
            planeControls.style.display = show ? 'flex' : 'none';
        }

        function setPlaneFilter(enabled) {
            state.planeFilterEnabled = !!enabled;
            if (planeToggleBtn) planeToggleBtn.classList.toggle('active', state.planeFilterEnabled);
            if (state.planeFilterEnabled) {
                // Prefer full plane range from img_dim.n_planes; fallback to computed range
                let minP = 0, maxP = 0;
                if (state.geom && state.geom.imgDim && typeof state.geom.imgDim.n_planes === 'number') {
                    const n = Math.max(0, parseInt(state.geom.imgDim.n_planes, 10) || 0);
                    minP = 0;
                    maxP = Math.max(0, n - 1);
                } else if (state.geom && state.geom.planeId && state.planeRange) {
                    minP = state.planeRange.min;
                    maxP = state.planeRange.max;
                }
                if (planeSlider) {
                    planeSlider.min = String(minP);
                    planeSlider.max = String(maxP);
                    if (state.selectedPlane === null) {
                        // Default to middle plane when enabling
                        state.selectedPlane = Math.floor((minP + maxP) / 2);
                    }
                    planeSlider.value = String(state.selectedPlane);
                }
                if (planeLabel) planeLabel.textContent = `Plane: ${state.selectedPlane}`;
            } else {
                state.selectedPlane = null;
            }
            updatePlaneFooterVisibility();
            window.pciSeq.render();
        }

        if (planeToggleBtn) {
            planeToggleBtn.addEventListener('click', () => {
                if (!state.geom.is3D || !state.geom.planeId) {
                    window.pciSeq.showUserNotice('Plane filter is available only for 3D data.');
                    return;
                }
                setPlaneFilter(!state.planeFilterEnabled);
            });
            // Always show the toggle; disable it when not 3D
            planeToggleBtn.disabled = !(state.geom.is3D);
        }
        if (planeSlider) {
            planeSlider.addEventListener('input', (e) => {
                state.selectedPlane = parseInt(e.target.value, 10);
                if (planeLabel) planeLabel.textContent = `Plane: ${state.selectedPlane}`;
                window.pciSeq.render();
            });
        }
    }

    // Setup threshold slider for changes view
    function setupThresholdSlider() {
        const thresholdSlider = document.getElementById('threshold-slider');
        const thresholdValueDisplay = document.getElementById('threshold-value');

        if (thresholdSlider && thresholdValueDisplay) {
            thresholdSlider.addEventListener('input', (e) => {
                state.changeThreshold = parseFloat(e.target.value);
                thresholdValueDisplay.textContent = state.changeThreshold.toFixed(1);

                // Recalculate changed cells with new threshold, but DON'T update previous prob
                if (state.cells.length > 0 && state.previousProb) {
                    window.pciSeq.detection.detectChangedCells(false);
                    window.pciSeq.render();
                }
            });
        }
    }

    // Setup legend filter with keyboard shortcut
    function setupLegendFilter() {
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
                window.pciSeq.updateLegend();
            });
            legendFilterInput.addEventListener('keydown', (e) => {
                if (e.key === 'Escape') {
                    legendFilterInput.value = '';
                    state.legendFilter = '';
                    window.pciSeq.updateLegend();
                }
            });
        }
    }

    // Setup window resize handler
    function setupWindowResize() {
        window.addEventListener('resize', () => {
            if (state.deckgl) {
                state.deckgl.setProps({
                    width: '100%',
                    height: '100%'
                });
            }
        });
    }

    // Export initialization function
    window.pciSeq.uiControls = {
        initializeControls: initializeControls
    };

})();