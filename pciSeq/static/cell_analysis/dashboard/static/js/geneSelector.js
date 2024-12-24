import { PLOT_CONFIG } from './plotConfig.js';

export class GeneSelector {
    constructor(containerId, geneNames, onChange) {
        // Core properties
        this.containerId = containerId;
        this.geneNames = geneNames;
        this.selectedGenes = new Set(geneNames); // Start with all genes selected
        this.onChange = onChange;
        
        // Search and filter properties
        this.searchTerm = '';
        this.container = d3.select(`#${this.containerId}`);
        
        // Initialize the component
        this.initializeSelector();
    }

    initializeSelector() {
        // Clear any existing content
        this.container.html('');
        
        // Create header section
        const header = this.container.append('div')
            .attr('class', 'gene-selector-header');
            
        header.append('h3')
            .text('Gene Panel');
            
        // Create control buttons
        const controls = header.append('div')
            .attr('class', 'gene-selector-controls');
            
        controls.append('button')
            .attr('class', 'control-button')
            .text('Select All')
            .on('click', () => this.selectAll());
            
        controls.append('button')
            .attr('class', 'control-button')
            .text('Deselect All')
            .on('click', () => this.deselectAll());

        // Add search input
        const searchContainer = this.container.append('div')
            .attr('class', 'gene-search-container');

        searchContainer.append('input')
            .attr('type', 'text')
            .attr('class', 'gene-search-input')
            .attr('placeholder', 'Search genes...')
            .on('input', (event) => {
                this.searchTerm = event.target.value.toLowerCase();
                this.updateVisibleCheckboxes();
            });

        // Create checkbox container
        this.checkboxContainer = this.container.append('div')
            .attr('class', 'gene-checkbox-grid');

        // Initialize checkboxes
        this.createCheckboxes();
    }

    createCheckboxes() {
        // Create checkbox for each gene
        const items = this.checkboxContainer.selectAll('div')
            .data(this.geneNames)
            .enter()
            .append('div')
            .attr('class', 'gene-checkbox-item')
            .style('display', d => 
                this.searchTerm ? 
                    (d.toLowerCase().includes(this.searchTerm) ? 'flex' : 'none') 
                    : 'flex'
            );

        // Add checkbox input
        items.append('input')
            .attr('type', 'checkbox')
            .attr('id', d => `gene-${this.containerId}-${d}`)
            .attr('checked', d => this.selectedGenes.has(d))
            .on('change', (event, d) => {
                if (event.target.checked) {
                    this.selectedGenes.add(d);
                } else {
                    this.selectedGenes.delete(d);
                }
                this.onChange(Array.from(this.selectedGenes));
            });

        // Add label
        items.append('label')
            .attr('for', d => `gene-${this.containerId}-${d}`)
            .attr('class', 'gene-label')
            .text(d => d);
    }

    updateVisibleCheckboxes() {
        this.checkboxContainer.selectAll('.gene-checkbox-item')
            .style('display', d => 
                this.searchTerm ? 
                    (d.toLowerCase().includes(this.searchTerm) ? 'flex' : 'none') 
                    : 'flex'
            );
    }

    selectAll() {
        this.selectedGenes = new Set(this.geneNames);
        this.updateCheckboxes(true);
        this.onChange(Array.from(this.selectedGenes));
    }

    deselectAll() {
        this.selectedGenes.clear();
        this.updateCheckboxes(false);
        this.onChange(Array.from(this.selectedGenes));
    }

    updateCheckboxes(checked) {
        this.checkboxContainer.selectAll('input[type="checkbox"]')
            .property('checked', checked);
    }

    // Method to programmatically select specific genes
    selectGenes(genes) {
        genes.forEach(gene => {
            if (this.geneNames.includes(gene)) {
                this.selectedGenes.add(gene);
            }
        });
        this.updateCheckboxes();
        this.onChange(Array.from(this.selectedGenes));
    }

    // Method to get currently selected genes
    getSelectedGenes() {
        return Array.from(this.selectedGenes);
    }

    // Method to handle errors
    handleError(error) {
        console.error('GeneSelector Error:', error);
        this.container.html(`
            <div class="error-message">
                An error occurred: ${error.message}
                <button onclick="location.reload()">Reload</button>
            </div>
        `);
    }
}
