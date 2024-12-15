export class GeneSelector {
    constructor(containerId, genes, onSelectionChange) {
        this.container = document.getElementById(containerId);
        this.genes = genes;
        this.onSelectionChange = onSelectionChange;
        this.selectedGenes = new Set(genes); // Start with all genes selected
        this.init();
    }

    init() {
        // Initialize select all checkbox - unchecked by default
        const selectAllCheckbox = this.container.querySelector('#select-all-genes');
        selectAllCheckbox.checked = false; // Start unchecked
        selectAllCheckbox.addEventListener('change', () => this.handleSelectAll(!selectAllCheckbox.checked));

        const checkboxContainer = this.container.querySelector('#gene-checkboxes');

        // Create checkbox for each gene
        this.genes.forEach(geneName => {
            const wrapper = document.createElement('div');
            wrapper.className = 'gene-checkbox-wrapper';

            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.id = `gene-${geneName}`;
            checkbox.checked = true; // Start with all genes checked
            checkbox.addEventListener('change', () => this.handleCheckboxChange(geneName));

            const label = document.createElement('label');
            label.htmlFor = `gene-${geneName}`;
            label.textContent = geneName;

            wrapper.appendChild(checkbox);
            wrapper.appendChild(label);
            checkboxContainer.appendChild(wrapper);
        });
    }

    handleSelectAll(shouldSelect) {
        // Update all checkboxes
        this.genes.forEach(geneName => {
            const checkbox = document.getElementById(`gene-${geneName}`);
            checkbox.checked = shouldSelect;
        });

        // Update selected genes set
        this.selectedGenes = new Set(shouldSelect ? this.genes : []);

        // Update toggle label - note the inverted logic
        const toggleLabel = this.container.querySelector('.toggle-label');
        toggleLabel.textContent = shouldSelect ? 'Unselect All' : 'Select All';

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

        // Update select all checkbox state - note the inverted logic
        const selectAllCheckbox = this.container.querySelector('#select-all-genes');
        const toggleLabel = this.container.querySelector('.toggle-label');

        if (this.selectedGenes.size === this.genes.length) {
            selectAllCheckbox.checked = false;
            toggleLabel.textContent = 'Unselect All';
        } else if (this.selectedGenes.size === 0) {
            selectAllCheckbox.checked = true;
            toggleLabel.textContent = 'Select All';
        }

        // Notify plots of the change
        this.onSelectionChange(Array.from(this.selectedGenes));
    }

    getSelectedGenes() {
        return Array.from(this.selectedGenes);
    }
}
