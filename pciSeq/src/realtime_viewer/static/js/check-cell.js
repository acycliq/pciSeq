/**
 * Check Cell Module
 * Handles cell diagnostics via Ctrl+Click interaction
 */

(function() {
    'use strict';

    window.pciSeq = window.pciSeq || {};
    const state = window.pciSeq.state;

    // Current selected cell for diagnostics
    let currentCellLabel = null;
    let currentCellClass = null;
    let currentComparisonClass = 'Zero';

    /**
     * Initialize check cell functionality
     */
    function initialize() {
        console.log('=== Initializing check cell module ===');

        // Set up close button
        const closeBtn = document.getElementById('drawer-close');
        if (closeBtn) {
            closeBtn.addEventListener('click', closeDrawer);
            console.log('Close button handler attached');
        } else {
            console.error('drawer-close button not found!');
        }

        // Set up comparison class dropdown
        const selectEl = document.getElementById('comparison-class-select');
        if (selectEl) {
            selectEl.addEventListener('change', (e) => {
                currentComparisonClass = e.target.value;
                if (currentCellLabel !== null) {
                    requestCheckCell(currentCellLabel, currentComparisonClass);
                }
            });
            console.log('Comparison dropdown handler attached');
        } else {
            console.error('comparison-class-select not found!');
        }

        console.log('Check cell module initialized successfully');
    }

    /**
     * Populate the comparison class dropdown with available classes
     */
    function populateClassDropdown() {
        const selectEl = document.getElementById('comparison-class-select');
        if (!selectEl || !state.cellClassNames) return;

        // Clear existing options
        selectEl.innerHTML = '';

        // state.cellClassNames is an object {index: name}, so iterate over values
        const classNames = Object.values(state.cellClassNames);

        // Add all class names
        classNames.forEach((className) => {
            const option = document.createElement('option');
            option.value = className;
            option.textContent = className;
            // Select the current comparison class, not just 'Zero'
            if (className === currentComparisonClass) {
                option.selected = true;
            }
            selectEl.appendChild(option);
        });

        // Also ensure the select element's value property is set, as `option.selected`
        // might not be enough in all browsers to reflect the state.
        selectEl.value = currentComparisonClass;

        console.log(`Populated dropdown with ${classNames.length} classes, selected: ${currentComparisonClass}`);
    }

    /**
     * Handle cell click event from deck.gl
     * @param {Object} info - Click info from deck.gl
     * @param {Object} event - Mouse event
     */
    function handleCellClick(info, event) {
        console.log('=== handleCellClick called ===', {
            hasInfo: !!info,
            hasObject: !!(info && info.object),
            ctrlKey: event.ctrlKey,
            metaKey: event.metaKey
        });

        // Only respond to Ctrl+Click
        if (!event.ctrlKey && !event.metaKey) {
            console.log('Not a Ctrl+Click, ignoring');
            return;
        }

        if (!info || !info.object) {
            console.log('No cell clicked');
            return;
        }

        const cell = info.object;
        console.log('Cell object:', cell);

        // Get cell label (ID) - assuming cell.id is the label
        // Note: may need to map this to the original cell label
        currentCellLabel = cell.id;
        currentCellClass = state.cellClassNames[cell.class] || `Class ${cell.class}`;

        console.log(`Ctrl+Click on cell ${currentCellLabel}, class: ${currentCellClass}`);

        // Make sure dropdown is populated
        populateClassDropdown();

        // Open drawer and request data
        openDrawer();
        requestCheckCell(currentCellLabel, currentComparisonClass);
    }

    /**
     * Open the left side drawer panel
     */
    function openDrawer() {
        const drawer = document.getElementById('bottom-drawer');
        if (drawer) {
            drawer.classList.add('open');

            // Update title
            const title = document.getElementById('drawer-title');
            if (title && currentCellLabel !== null) {
                title.innerHTML = `Cell ${currentCellLabel}:<br/>${currentCellClass}`;
            }
        }
    }

    /**
     * Close the bottom drawer panel
     */
    function closeDrawer() {
        const drawer = document.getElementById('bottom-drawer');
        if (drawer) {
            drawer.classList.remove('open');
        }
        currentCellLabel = null;
        currentCellClass = null;
    }

    /**
     * Request check_cell data from server via WebSocket
     * @param {number} cellLabel - Cell ID/label
     * @param {string} comparisonClass - Class name to compare against
     */
    function requestCheckCell(cellLabel, comparisonClass) {
        // Get socket from global namespace (set by socket-handlers.js)
        const socket = window.pciSeq.socket || state.socket;

        if (!socket) {
            console.error('Cannot request check_cell: socket not available');
            showErrorInCharts('Socket not available');
            return;
        }

        if (!socket.connected) {
            console.error('Cannot request check_cell: socket not connected');
            showErrorInCharts('Socket not connected. Please wait for connection.');
            return;
        }

        console.log(`Requesting check_cell for cell ${cellLabel} vs ${comparisonClass}`);

        // Show loading state in both chart areas
        showLoadingInCharts();

        // Send request to server
        socket.emit('request_check_cell', {
            cell_label: cellLabel,
            comparison_class: comparisonClass
        });
    }

    /**
     * Show loading indicator in chart areas
     */
    function showLoadingInCharts() {
        const leftDiv = document.getElementById('check-cell-chart-left');
        const rightDiv = document.getElementById('check-cell-chart-right');
        const loadingHTML = '<div style="display: flex; align-items: center; justify-content: center; height: 100%; color: var(--muted); font-size: 12px;">Calculating...</div>';

        if (leftDiv) leftDiv.innerHTML = loadingHTML;
        if (rightDiv) rightDiv.innerHTML = loadingHTML;
    }

    /**
     * Show error message in chart areas
     */
    function showErrorInCharts(message) {
        const leftDiv = document.getElementById('check-cell-chart-left');
        const rightDiv = document.getElementById('check-cell-chart-right');
        const errorHTML = `<div style="display: flex; align-items: center; justify-content: center; height: 100%; color: var(--danger); font-size: 11px; text-align: center; padding: 20px;">${message}</div>`;

        if (leftDiv) leftDiv.innerHTML = errorHTML;
        if (rightDiv) rightDiv.innerHTML = errorHTML;
    }

    /**
     * Handle check_cell response from server
     * @param {Object} data - Response data containing gene comparison data
     */
    function handleCheckCellResponse(data) {
        console.log('Received check_cell response:', data);

        if (data.error) {
            showErrorInCharts(`Error: ${data.error}`);
            return;
        }

        // Render the D3 chart
        renderCheckCellChart(data);
    }

    /**
     * Render the gene comparison chart using D3
     * @param {Object} data - Chart data from server
     */
    function renderCheckCellChart(data) {
        const leftDiv = document.getElementById('check-cell-chart-left');
        const rightDiv = document.getElementById('check-cell-chart-right');
        if (!leftDiv || !rightDiv) return;

        // Clear previous charts
        leftDiv.innerHTML = '';
        rightDiv.innerHTML = '';

        // Extract data
        const topGenes = data.top_genes || [];
        const bottomGenes = data.bottom_genes || [];
        const pciSeqClass = data.pciseq_class || '';
        const userClass = data.user_class || '';
        const cellLabel = data.cell_label || '';
        const topSum = data.top_sum || 0;
        const bottomSum = data.bottom_sum || 0;
        const geneData = data.gene_expression_data || null;

        // Setup dimensions for each chart - tight margins to maximize chart area
        // Reduce empty space above the title; keep bottom for x-axis labels
        const margin = {top: 10, right: 10, bottom: 60, left: 45};
        const containerWidth = leftDiv.clientWidth;
        const containerHeight = leftDiv.clientHeight;
        const chartWidth = containerWidth - margin.left - margin.right;
        const chartHeight = containerHeight - margin.top - margin.bottom;

        // Color scheme
        const colorTop = '#87CEEB'; // skyblue
        const colorBottom = '#FA8072'; // salmon

        // Render left chart (top genes for pciSeq class) - positive bars
        renderSingleBarChart(
            leftDiv,
            topGenes,
            margin,
            chartWidth,
            chartHeight,
            colorTop,
            `Cell ${cellLabel} - Top 10 contr for class:\n${pciSeqClass} (Sum: ${topSum.toFixed(2)})`,
            false
        );

        // Render right chart (bottom genes for user class) - negative bars
        renderSingleBarChart(
            rightDiv,
            bottomGenes,
            margin,
            chartWidth,
            chartHeight,
            colorBottom,
            `Cell ${cellLabel} - Top 10 contr for class:\n${userClass} (Sum: ${bottomSum.toFixed(2)})`,
            true
        );

        // Render gene expression data table if available
        if (geneData && geneData.length > 0) {
            renderGeneTable(geneData, pciSeqClass, userClass, topGenes, bottomGenes);
        }
    }

    /**
     * Render a single bar chart in its container
     */
    function renderSingleBarChart(container, genes, margin, chartWidth, chartHeight, color, title, negative) {
        // Create SVG
        const svg = d3.select(container)
            .append('svg')
            .attr('width', chartWidth + margin.left + margin.right)
            .attr('height', chartHeight + margin.top + margin.bottom);

        const chartGroup = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        // Reserve a small gap between the multi-line title and the plot area
        const titleGap = 20; // px

        // X scale (band scale for gene names)
        const xScale = d3.scaleBand()
            .domain(d3.range(genes.length))
            .range([0, chartWidth])
            .padding(0.2);

        // Y scale (linear scale for values) over the available plot height (excludes titleGap)
        let yDomain;
        if (negative) {
            const minVal = d3.min(genes, d => d.value);
            const minNeg = (minVal !== undefined && minVal !== null) ? Math.min(0, minVal) : -1;
            yDomain = [minNeg * 1.1, 0];
        } else {
            const maxVal = d3.max(genes, d => d.value);
            const maxPos = (maxVal !== undefined && maxVal !== null) ? Math.max(0, maxVal) : 1;
            yDomain = [0, maxPos * 1.1];
        }
        const plotHeight = Math.max(1, chartHeight - titleGap);
        const yScale = d3.scaleLinear()
            .domain(yDomain)
            .range([plotHeight, 0]);

        // Title (supports multi-line with \n)
        const titleLines = title.split('\n');
        const titleGroup = chartGroup.append('text')
            .attr('x', chartWidth / 2)
            .attr('y', 0) // keep title within the SVG group
            .attr('text-anchor', 'middle')
            .style('font-size', '10px')
            .style('font-weight', '600')
            .style('fill', 'var(--text)');

        titleLines.forEach((line, i) => {
            titleGroup.append('tspan')
                .attr('x', chartWidth / 2)
                .attr('dy', i === 0 ? 0 : '1.2em')
                .text(line);
        });

        // Render bars and axes inside a plot group offset by titleGap
        const plotGroup = chartGroup.append('g').attr('transform', `translate(0, ${titleGap})`);
        renderBarsAndAxes(plotGroup, genes, xScale, yScale, color, plotHeight, negative);
    }

    /**
     * Helper to render bars and axes
     */
    function renderBarsAndAxes(chartGroup, genes, xScale, yScale, color, chartHeight, negative) {

        // Bars relative to zero baseline for negative charts
        const zeroY = yScale(0);
        chartGroup.selectAll('.bar')
            .data(genes)
            .join('rect')
            .attr('class', 'bar')
            .attr('x', (d, i) => xScale(i))
            .attr('y', d => negative ? zeroY : yScale(d.value))
            .attr('width', xScale.bandwidth())
            .attr('height', d => {
                if (negative) {
                    const vy = yScale(d.value);
                    return Math.max(0, vy - zeroY);
                } else {
                    return Math.max(0, zeroY - yScale(d.value));
                }
            })
            .attr('fill', color)
            .attr('opacity', 0.8)
            .on('mouseover', function(event, d) {
                d3.select(this).attr('opacity', 1);
                showTooltip(event, d);
            })
            .on('mouseout', function() {
                d3.select(this).attr('opacity', 0.8);
                hideTooltip();
            });

        // Zero baseline
        chartGroup.append('line')
            .attr('x1', 0)
            .attr('x2', xScale.range()[1])
            .attr('y1', zeroY)
            .attr('y2', zeroY)
            .attr('stroke', 'var(--border)')
            .attr('stroke-width', 1);

        // X-axis with gene names at the bottom
        const xAxis = d3.axisBottom(xScale)
            .tickFormat((_, i) => (genes[i] && genes[i].gene) ? genes[i].gene : '');

        const xAxisGroup = chartGroup.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0, ${yScale.range()[0]})`)
            .call(xAxis);

        // Style and rotate tick labels for readability
        xAxisGroup.selectAll('text')
            .style('text-anchor', 'end')
            .style('fill', 'var(--text)')
            .style('font-size', '10px')
            .attr('dx', '-0.5em')
            .attr('dy', '0.15em')
            .attr('transform', 'rotate(-45)');

        // Y-axis
        const yAxis = d3.axisLeft(yScale).ticks(5);
        chartGroup.append('g')
            .attr('class', 'y-axis')
            .call(yAxis)
            .style('color', 'var(--muted)');

        // Y-axis label
        chartGroup.append('text')
            .attr('transform', 'rotate(-90)')
            .attr('x', -chartHeight / 2)
            .attr('y', -32)
            .attr('text-anchor', 'middle')
            .style('font-size', '10px')
            .style('fill', 'var(--muted)')
            .text('Log-Likelihood Difference');
    }

    /**
     * Show tooltip
     */
    function showTooltip(event, d) {
        const tooltip = d3.select('body')
            .append('div')
            .attr('class', 'chart-tooltip visible')
            .style('position', 'fixed')
            .style('z-index', '3001')
            .style('pointer-events', 'none')
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 10) + 'px')
            .html(`
                <strong>${d.gene}</strong><br/>
                Contribution: ${d.value.toFixed(3)}
            `);
    }

    /**
     * Hide tooltip
     */
    function hideTooltip() {
        d3.selectAll('.chart-tooltip').remove();
    }

    /**
     * Render gene expression data table
     * @param {Array} geneData - Array of gene expression data objects
     * @param {string} pciSeqClass - pciSeq assigned class name
     * @param {string} userClass - User comparison class name
     */
    function renderGeneTable(geneData, pciSeqClass, userClass, topGenes, bottomGenes) {
        const tableEl = document.getElementById('check-cell-table');
        if (!tableEl) return;

        // Clear existing table
        tableEl.innerHTML = '';

        // Create table header
        const thead = document.createElement('thead');
        const headerRow = document.createElement('tr');

        const headers = [
            'Gene',
            `Mean counts (cells typed as ${pciSeqClass})`,
            `Mean counts (cells typed as ${userClass})`,
            'This cell: counts'
        ];

        headers.forEach(headerText => {
            const th = document.createElement('th');
            th.textContent = headerText;
            headerRow.appendChild(th);
        });

        thead.appendChild(headerRow);
        tableEl.appendChild(thead);

        // Create table body
        const tbody = document.createElement('tbody');

        // Index geneData by gene for quick lookups
        const byGene = new Map(geneData.map(r => [r.gene, r]));

        // Helper to append a subheader row
        function appendSubheader(text) {
            const sr = document.createElement('tr');
            sr.className = 'subheader';
            const td = document.createElement('td');
            td.colSpan = 4;
            td.textContent = text;
            sr.appendChild(td);
            tbody.appendChild(sr);
        }

        // Render top 10 (supports assigned)
        if (Array.isArray(topGenes) && topGenes.length > 0) {
            appendSubheader(`Top 10 genes supporting Assigned (${pciSeqClass})`);
            topGenes.forEach(g => {
                const row = byGene.get(g.gene);
                if (!row) return;
                const tr = document.createElement('tr');

                // Gene name
                const tdGene = document.createElement('td');
                tdGene.className = 'gene-name assigned';
                tdGene.textContent = row.gene;
                tr.appendChild(tdGene);

                // Means
                const tdP = document.createElement('td');
                tdP.textContent = Number(row.mean_expr_pciseq).toFixed(3);
                tr.appendChild(tdP);
                const tdU = document.createElement('td');
                tdU.textContent = Number(row.mean_expr_user).toFixed(3);
                tr.appendChild(tdU);

                // This cell count
                const tdC = document.createElement('td');
                tdC.textContent = Number(row.gene_count).toFixed(3);
                tr.appendChild(tdC);

                tbody.appendChild(tr);
            });
        }

        // Render bottom 10 (supports comparison)
        if (Array.isArray(bottomGenes) && bottomGenes.length > 0) {
            appendSubheader(`Top 10 genes supporting Comparison (${userClass})`);
            bottomGenes.forEach(g => {
                const row = byGene.get(g.gene);
                if (!row) return;
                const tr = document.createElement('tr');

                // Gene name
                const tdGene = document.createElement('td');
                tdGene.className = 'gene-name comparison';
                tdGene.textContent = row.gene;
                tr.appendChild(tdGene);

                // Means
                const tdP = document.createElement('td');
                tdP.textContent = Number(row.mean_expr_pciseq).toFixed(3);
                tr.appendChild(tdP);
                const tdU = document.createElement('td');
                tdU.textContent = Number(row.mean_expr_user).toFixed(3);
                tr.appendChild(tdU);

                // This cell count
                const tdC = document.createElement('td');
                tdC.textContent = Number(row.gene_count).toFixed(3);
                tr.appendChild(tdC);

                tbody.appendChild(tr);
            });
        }

        tableEl.appendChild(tbody);

        console.log(`Rendered gene table with ${geneData.length} genes`);
    }

    // Export public API
    window.pciSeq.checkCell = {
        initialize,
        handleCellClick,
        handleCheckCellResponse,
        populateClassDropdown
    };

})();
