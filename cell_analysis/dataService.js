// Save this as: static/js/api/dataService.js

export class DataService {
    static async fetchCellData(cellId) {
        try {
            const response = await fetch(`/api/cell-data/${cellId}`);
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            const data = await response.json();
            return this.validateData(data);
        } catch (error) {
            console.error('Error fetching cell data:', error);
            throw error;
        }
    }

    static validateData(data) {
        // Validate required data properties
        const requiredProps = [
            'cell_num',
            'class_names',
            'gene_names',
            'user_class',
            'assigned_class',
            'class_probs',
            'contr'
        ];

        for (const prop of requiredProps) {
            if (!(prop in data)) {
                throw new Error(`Missing required property: ${prop}`);
            }
        }

        // Validate data structure
        if (!Array.isArray(data.class_names) || !Array.isArray(data.gene_names)) {
            throw new Error('class_names and gene_names must be arrays');
        }

        if (typeof data.contr !== 'object') {
            throw new Error('contr must be an object');
        }

        return data;
    }
}
