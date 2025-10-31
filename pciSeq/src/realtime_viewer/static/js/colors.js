/**
 * Color Management Module
 * Handles color palette generation, custom color schemes, and color conversions
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Generate distinct colors for up to 65 classes
    function generateColorPalette() {
        const numColors = 65;
        const hueStep = 360 / numColors;

        for (let i = 0; i < numColors; i++) {
            const hue = (i * hueStep) % 360;
            const saturation = 70 + (i % 3) * 10; // Vary saturation slightly
            const lightness = 50 + (i % 2) * 10;  // Vary lightness slightly
            state.cellClassColors[i] = hslToRgb(hue, saturation, lightness);
        }
    }

    // HSL to RGB conversion
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

    // Convert hex color to RGB using D3.js
    function hexToRgb(hex) {
        const color = d3.rgb(hex);
        return [color.r, color.g, color.b];
    }

    // Apply custom color scheme to cell classes
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

    // Load custom colors from user-provided scheme
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
            window.pciSeq.updateLegend();
            window.pciSeq.render();
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

    // Handle color file upload
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

    // Export functions
    window.pciSeq.colors = {
        generateColorPalette: generateColorPalette,
        applyColorScheme: applyColorScheme,
        loadCustomColors: loadCustomColors,
        handleColorFileUpload: handleColorFileUpload
    };

})();