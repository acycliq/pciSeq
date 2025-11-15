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
            // Select Zero by default
            if (className === 'Zero') {
                option.selected = true;
            }
            selectEl.appendChild(option);
        });

        console.log(`Populated dropdown with ${classNames.length} classes`);
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
     * Open the bottom drawer panel
     */
    function openDrawer() {
        const drawer = document.getElementById('bottom-drawer');
        if (drawer) {
            drawer.classList.add('open');

            // Update title
            const title = document.getElementById('drawer-title');
            if (title && currentCellLabel !== null) {
                title.textContent = `Cell ${currentCellLabel} (${currentCellClass})`;
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
            const chartDiv = document.getElementById('check-cell-chart');
            if (chartDiv) {
                chartDiv.innerHTML = '<div style="text-align: center; padding: 40px; color: var(--danger);">Socket not available</div>';
            }
            return;
        }

        if (!socket.connected) {
            console.error('Cannot request check_cell: socket not connected');
            const chartDiv = document.getElementById('check-cell-chart');
            if (chartDiv) {
                chartDiv.innerHTML = '<div style="text-align: center; padding: 40px; color: var(--danger);">Socket not connected. Please wait for connection.</div>';
            }
            return;
        }

        console.log(`Requesting check_cell for cell ${cellLabel} vs ${comparisonClass}`);

        // Show loading state in chart area
        const chartDiv = document.getElementById('check-cell-chart');
        if (chartDiv) {
            chartDiv.innerHTML = '<div style="text-align: center; padding: 40px; color: var(--muted);">Loading...</div>';
        }

        // Send request to server
        socket.emit('request_check_cell', {
            cell_label: cellLabel,
            comparison_class: comparisonClass
        });
    }

    /**
     * Handle check_cell response from server
     * @param {Object} data - Response data containing gene comparison data
     */
    function handleCheckCellResponse(data) {
        console.log('Received check_cell response:', data);

        if (data.error) {
            const chartDiv = document.getElementById('check-cell-chart');
            if (chartDiv) {
                chartDiv.innerHTML = `<div style="text-align: center; padding: 40px; color: var(--danger);">Error: ${data.error}</div>`;
            }
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
        const chartDiv = document.getElementById('check-cell-chart');
        if (!chartDiv) return;

        // Clear previous chart
        chartDiv.innerHTML = '';

        // Extract data
        const topGenes = data.top_genes || [];
        const bottomGenes = data.bottom_genes || [];
        const pciSeqClass = data.pciseq_class || '';
        const userClass = data.user_class || '';
        const topSum = data.top_sum || 0;
        const bottomSum = data.bottom_sum || 0;
        const geneData = data.gene_expression_data || null;

        // Setup dimensions
        const margin = {top: 40, right: 20, bottom: 80, left: 60};
        const containerWidth = chartDiv.clientWidth;
        const containerHeight = chartDiv.clientHeight;
        const chartWidth = Math.floor((containerWidth - margin.left - margin.right) / 2) - 10;
        const chartHeight = containerHeight - margin.top - margin.bottom;

        // Create SVG
        const svg = d3.select(chartDiv)
            .append('svg')
            .attr('width', containerWidth)
            .attr('height', containerHeight);

        // Create two chart groups (side by side)
        const leftChart = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);

        const rightChart = svg.append('g')
            .attr('transform', `translate(${margin.left + chartWidth + 40},${margin.top})`);

        // X scales (band scale for gene names - VERTICAL bars)
        const xScaleTop = d3.scaleBand()
            .domain(d3.range(topGenes.length))
            .range([0, chartWidth])
            .padding(0.2);

        const xScaleBottom = d3.scaleBand()
            .domain(d3.range(bottomGenes.length))
            .range([0, chartWidth])
            .padding(0.2);

        // Y scales (linear scale for values - VERTICAL bars)
        const topMax = d3.max(topGenes, d => d.value) || 1;
        const bottomMax = d3.max(bottomGenes, d => Math.abs(d.value)) || 1;

        const yScaleTop = d3.scaleLinear()
            .domain([0, topMax * 1.1])
            .range([chartHeight, 0]);

        const yScaleBottom = d3.scaleLinear()
            .domain([0, bottomMax * 1.1])
            .range([chartHeight, 0]);

        // Color scheme - match matplotlib skyblue and salmon
        const colorTop = '#87CEEB'; // skyblue
        const colorBottom = '#FA8072'; // salmon

        // Render top genes (left chart)
        renderBarChart(leftChart, topGenes, xScaleTop, yScaleTop, colorTop,
            `Top genes for ${pciSeqClass} (Sum: ${topSum.toFixed(2)})`);

        // Render bottom genes (right chart)
        renderBarChart(rightChart, bottomGenes, xScaleBottom, yScaleBottom, colorBottom,
            `Top genes for ${userClass} (Sum: ${Math.abs(bottomSum).toFixed(2)})`);

        // Render gene expression data table if available
        if (geneData && geneData.length > 0) {
            renderGeneTable(geneData, pciSeqClass, userClass);
        }
    }

    /**
     * Helper to render a single VERTICAL bar chart (like matplotlib)
     */
    function renderBarChart(chartGroup, genes, xScale, yScale, color, title, labelAlign) {
        // Title
        chartGroup.append('text')
            .attr('x', xScale.range()[1] / 2)
            .attr('y', -10)
            .attr('text-anchor', 'middle')
            .style('font-size', '13px')
            .style('font-weight', '600')
            .style('fill', 'var(--text)')
            .text(title);

        // Bars (VERTICAL - like matplotlib)
        // For bottom genes, show absolute values as positive bars
        chartGroup.selectAll('.bar')
            .data(genes)
            .join('rect')
            .attr('class', 'bar')
            .attr('x', (d, i) => xScale(i))
            .attr('y', d => yScale(Math.abs(d.value)))
            .attr('width', xScale.bandwidth())
            .attr('height', d => yScale.range()[0] - yScale(Math.abs(d.value)))
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

        // Gene labels on X-axis (bottom, rotated like matplotlib)
        chartGroup.selectAll('.label')
            .data(genes)
            .join('text')
            .attr('class', 'label')
            .attr('x', (d, i) => xScale(i) + xScale.bandwidth() / 2)
            .attr('y', yScale.range()[1] + 10)
            .attr('dy', '0.35em')
            .attr('text-anchor', 'end')
            .attr('transform', (d, i) => `rotate(-45, ${xScale(i) + xScale.bandwidth() / 2}, ${yScale.range()[1] + 10})`)
            .style('font-size', '10px')
            .style('fill', 'var(--text)')
            .text(d => d.gene);

        // Y-axis
        const yAxis = d3.axisLeft(yScale).ticks(5);
        chartGroup.append('g')
            .attr('class', 'y-axis')
            .call(yAxis)
            .style('color', 'var(--muted)');

        // Y-axis label
        chartGroup.append('text')
            .attr('transform', 'rotate(-90)')
            .attr('x', -yScale.range()[1] / 2)
            .attr('y', -50)
            .attr('text-anchor', 'middle')
            .style('font-size', '11px')
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
    function renderGeneTable(geneData, pciSeqClass, userClass) {
        const tableEl = document.getElementById('check-cell-table');
        if (!tableEl) return;

        // Clear existing table
        tableEl.innerHTML = '';

        // Create table header
        const thead = document.createElement('thead');
        const headerRow = document.createElement('tr');

        const headers = [
            'Gene',
            `Mean Expr (${pciSeqClass})`,
            `Mean Expr (${userClass})`,
            'Gene Count'
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

        geneData.forEach(row => {
            const tr = document.createElement('tr');

            // Gene name
            const tdGene = document.createElement('td');
            tdGene.className = 'gene-name';
            tdGene.textContent = row.gene;
            tr.appendChild(tdGene);

            // Mean expression for pciSeq class
            const tdPciSeq = document.createElement('td');
            tdPciSeq.textContent = row.mean_expr_pciseq.toFixed(3);
            tr.appendChild(tdPciSeq);

            // Mean expression for user class
            const tdUser = document.createElement('td');
            tdUser.textContent = row.mean_expr_user.toFixed(3);
            tr.appendChild(tdUser);

            // Gene count
            const tdCount = document.createElement('td');
            tdCount.textContent = row.gene_count;
            tr.appendChild(tdCount);

            tbody.appendChild(tr);
        });

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