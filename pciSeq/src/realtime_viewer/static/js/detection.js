/**
 * Change Detection Module
 * Handles detection of cell probability changes and class changes with fade-out animations
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Detect cells that changed probability (for changes view)
    function detectChangedCells(updatePrevious = true) {
        // Clear changed cells set
        state.changedCells.clear();

        // Safety guard: if the number of cells changed between steps,
        // pause change detection for this step and reset the baseline.
        if (state.previousProb && state.previousProb.length !== state.cells.length) {
            const prevN = state.previousProb.length;
            const currN = state.cells.length;
            const msg = `Paused the "Changes" view for this step because the number of cells changed (${prevN} → ${currN}). This can happen after reconnecting or restarting. The view will resume automatically on the next step.`;
            console.warn(msg);
            window.pciSeq.showUserNotice(msg);
            // Reset baseline to current values and skip diffing for this step
            state.previousProb = state.cells.map(cell => cell.prob);
            window.pciSeq.updateChangesCount();
            return;
        }

        // On first iteration (no previous data to compare), all cells are "changed"
        if (!state.previousProb) {
            for (let i = 0; i < state.cells.length; i++) {
                state.changedCells.add(state.cells[i].id);
            }
            console.log(`First iteration (no previous data): marked all ${state.cells.length} cells as changed`);
        } else {
            // Compare with previous iteration
            const threshold = state.changeThreshold / 100; // Convert percentage to decimal

            for (let i = 0; i < state.cells.length; i++) {
                const currentProb = state.cells[i].prob;
                const prevProb = state.previousProb[i];

                // Calculate absolute change in probability
                const change = Math.abs(currentProb - prevProb);

                // Mark as changed if exceeds threshold
                if (change >= threshold) {
                    state.changedCells.add(state.cells[i].id);
                }
            }

            console.log(`Detected ${state.changedCells.size}/${state.cells.length} changed cells (threshold: ${state.changeThreshold}%)`);
        }

        // Only store current probability when new iteration arrives, not when threshold changes
        if (updatePrevious) {
            state.previousProb = state.cells.map(cell => cell.prob);
        }

        // Update changes count display
        window.pciSeq.updateChangesCount();
    }

    // Detect cells that changed CLASS (not just probability) for highlighting with fade-out
    function detectClassChanges() {
        // Safety: if no previous class data or cell count mismatch, skip detection
        if (!state.previousClass || state.previousClass.length !== state.cells.length) {
            if (state.previousClass) {
                console.log(`Skipping class change detection: cell count mismatch (${state.previousClass.length} → ${state.cells.length})`);
            }
            return;
        }

        const now = Date.now();
        let changedCount = 0;

        // Compare current class with previous class
        for (let i = 0; i < state.cells.length; i++) {
            const currentClass = state.cells[i].class;
            const prevClass = state.previousClass[i];

            // If class changed, add to highlight map with current timestamp
            if (currentClass !== prevClass) {
                state.classChangedCells.set(state.cells[i].id, now);
                changedCount++;
            }
        }

        console.log(`Detected ${changedCount} cells that changed class (iteration ${state.iteration})`);

        // Start animation loop if there are cells to animate and loop isn't already running
        if (state.classChangedCells.size > 0 && !state.animationFrameId) {
            startHighlightAnimation();
        }
    }

    // Animation loop for fade-out highlight effect
    function startHighlightAnimation() {
        function animationLoop() {
            const now = Date.now();
            let activeAnimations = 0;

            // Clean up expired animations
            for (const [cellId, timestamp] of state.classChangedCells.entries()) {
                const elapsed = now - timestamp;
                if (elapsed >= state.highlightFadeDuration) {
                    state.classChangedCells.delete(cellId);
                } else {
                    activeAnimations++;
                }
            }

            // Re-render to update radius scaling
            if (state.deckgl && state.cells.length > 0) {
                window.pciSeq.render();
            }

            // Continue animation if there are still active animations
            if (activeAnimations > 0) {
                state.animationFrameId = requestAnimationFrame(animationLoop);
            } else {
                // Stop animation loop
                state.animationFrameId = null;
                console.log('Highlight animation complete');
            }
        }

        // Start the loop
        state.animationFrameId = requestAnimationFrame(animationLoop);
        console.log(`Starting highlight animation for ${state.classChangedCells.size} cells`);
    }

    // Export functions
    window.pciSeq.detection = {
        detectChangedCells: detectChangedCells,
        detectClassChanges: detectClassChanges,
        startHighlightAnimation: startHighlightAnimation
    };

})();