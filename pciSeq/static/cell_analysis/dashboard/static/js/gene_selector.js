export class GeneSelector {
    constructor(containerId, genes, onSelectionChange) {
        this.container = document.getElementById(containerId);
        this.genes = genes;
        this.onSelectionChange = onSelectionChange;
        this.selectedGenes = new Set(genes); // Start with all genes selected
        this.init();
    }

    init() {
        const checkboxContainer = this.container.querySelector('#gene-checkboxes');

        // Create checkbox for each gene
        this.genes.forEach(geneName => {
            const wrapper = document.createElement('div');
            wrapper.className = 'gene-checkbox-wrapper';

            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.id = `gene-${geneName}`;
            checkbox.checked = true;
            checkbox.addEventListener('change', () => this.handleCheckboxChange(geneName));

            const label = document.createElement('label');
            label.htmlFor = `gene-${geneName}`;
            label.textContent = geneName;

            wrapper.appendChild(checkbox);
            wrapper.appendChild(label);
            checkboxContainer.appendChild(wrapper);
        });
    }

    handleCheckboxChange(geneName) {
        const checkbox = document.getElementById(`gene-${geneName}`);

        if (checkbox.checked) {
            this.selectedGenes.add(geneName);
        } else {
            this.selectedGenes.delete(geneName);
        }

        // Notify plots of the change
        this.onSelectionChange(Array.from(this.selectedGenes));
    }

    getSelectedGenes() {
        return Array.from(this.selectedGenes);
    }
}
