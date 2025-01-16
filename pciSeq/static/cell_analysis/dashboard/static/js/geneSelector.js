import { PLOT_CONFIG } from './plotConfig.js';

export class GeneSelector {
    constructor(containerId, genes, onSelectionChange) {
        this.container = document.getElementById(containerId);
        this.genes = genes;
        this.selectedGenes = new Set(genes);
        this.onSelectionChange = onSelectionChange;
        this.filteredGenes = [...genes];
        this.searchTerm = '';
        this.render();
    }

    render() {
        this.container.innerHTML = `
            <div class="gene-selector-header">
                <h3>Gene Selection</h3>
                <div class="select-controls">
                    <label>
                        <span class="toggle-label">Unselect All</span>
                        <input type="checkbox" id="select-all-genes" checked>
                    </label>
                </div>
            </div>
            <div class="search-box">
                <input type="text" 
                       id="gene-search" 
                       placeholder="Search genes..."
                       value="${this.searchTerm}"
                       autocomplete="off">
            </div>
            <div class="gene-list">
                ${this.filteredGenes.map(gene => `
                    <label>
                        <input type="checkbox" 
                               id="gene-${gene}" 
                               value="${gene}" 
                               ${this.selectedGenes.has(gene) ? 'checked' : ''}>
                        ${gene}
                    </label>
                `).join('')}
            </div>
        `;

        // Add event listeners
        const selectAllCheckbox = this.container.querySelector('#select-all-genes');
        selectAllCheckbox.addEventListener('change', (e) => {
            this.handleSelectAll(e.target.checked);
        });

        const searchInput = this.container.querySelector('#gene-search');
        searchInput.addEventListener('input', (e) => {
            this.searchTerm = e.target.value;
            this.handleSearch(this.searchTerm);
        });

        this.filteredGenes.forEach(gene => {
            const checkbox = document.getElementById(`gene-${gene}`);
            checkbox.addEventListener('change', () => this.handleCheckboxChange(gene));
        });
    }

    handleSearch(searchTerm) {
        this.searchTerm = searchTerm;
        searchTerm = searchTerm.toLowerCase();
        this.filteredGenes = this.genes.filter(gene => 
            gene.toLowerCase().includes(searchTerm)
        );
        this.render();
        
        // Restore checkbox states after re-render
        this.selectedGenes.forEach(gene => {
            const checkbox = document.getElementById(`gene-${gene}`);
            if (checkbox) checkbox.checked = true;
        });

        // Restore search input focus
        const searchInput = this.container.querySelector('#gene-search');
        searchInput.focus();
        // Place cursor at the end of input
        searchInput.setSelectionRange(searchTerm.length, searchTerm.length);
    }

    handleSelectAll(checked) {
        // Update all checkboxes
        this.genes.forEach(geneName => {
            const checkbox = document.getElementById(`gene-${geneName}`);
            if (checkbox) checkbox.checked = checked;
        });

        // Update selected genes set
        this.selectedGenes = new Set(checked ? this.genes : []);

        // Update toggle label
        const toggleLabel = this.container.querySelector('.toggle-label');
        if (toggleLabel) {
            toggleLabel.textContent = checked ? 'Unselect All' : 'Select All';
        }

        // Notify plots of the change
        this.onSelectionChange(Array.from(this.selectedGenes));
    }

    handleCheckboxChange(geneName) {
        const checkbox = document.getElementById(`gene-${geneName}`);
        
        if (checkbox.checked) {
            this.selectedGenes.add(geneName);
        } else {
            this.selectedGenes.delete(geneName);
        }

        // Update select all checkbox state
        const selectAllCheckbox = this.container.querySelector('#select-all-genes');
        const toggleLabel = this.container.querySelector('.toggle-label');
        
        if (this.selectedGenes.size === this.genes.length) {
            selectAllCheckbox.checked = true;
            toggleLabel.textContent = 'Unselect All';
        } else {
            selectAllCheckbox.checked = false;
            toggleLabel.textContent = 'Select All';
        }

        // Notify plots of the change
        this.onSelectionChange(Array.from(this.selectedGenes));
    }

    // Method to programmatically select specific genes
    selectGenes(genes) {
        genes.forEach(gene => {
            if (this.genes.includes(gene)) {
                this.selectedGenes.add(gene);
            }
        });
        this.updateCheckboxes();
        this.onSelectionChange(Array.from(this.selectedGenes));
    }

    // Method to get currently selected genes
    getSelectedGenes() {
        return Array.from(this.selectedGenes);
    }

    // Method to handle errors
    handleError(error) {
        console.error('GeneSelector Error:', error);
        this.container.innerHTML = `
            <div class="error-message">
                An error occurred: ${error.message}
                <button onclick="location.reload()">Reload</button>
            </div>
        `;
    }
}
