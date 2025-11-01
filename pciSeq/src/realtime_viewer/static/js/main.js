/**
 * Main Initialization Module
 * Coordinates initialization of all modules and starts the application
 */

(function() {
    'use strict';

    // Initialize application when page loads
    window.addEventListener('load', function() {
        console.log('=== pciSeq Real-Time Viewer Initializing ===');

        // Initialize setup form (shown first)
        if (window.pciSeq.setup) {
            window.pciSeq.setup.initializeSetupForm();
        }

        // Initialize color palette
        window.pciSeq.colors.generateColorPalette();

        // Initialize deck.gl rendering (but visualization is hidden initially)
        window.pciSeq.rendering.initializeDeck();

        // Initialize UI controls
        window.pciSeq.uiControls.initializeControls();

        console.log('=== pciSeq Real-Time Viewer Ready ===');
    });

})();