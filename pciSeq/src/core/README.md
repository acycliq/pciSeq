# pciSeq Core Module

This module contains the core functionality of pciSeq, implementing the Variational Bayes algorithm for spatial transcriptomics analysis.

## Key Components

### 1. Main Algorithm (`main.py`)
- Implements the VarBayes class for iterative optimization
- Handles convergence and parameter updates
- Manages the overall analysis pipeline

### 2. Data Types (`datatypes.py`)
- Defines fundamental data structures:
  * Cells: Cell segmentation data
  * Genes: Gene expression data
  * Spots: RNA spot detection
  * SingleCell: scRNA-seq reference
  * CellType: Cell type classification

### 3. Summary Functions (`summary.py`)
- Processes analysis results
- Formats data for visualization
- Generates cell and spot summaries

### 4. Utilities (`utils.py`)
- Helper functions for computation
- File handling utilities
- Convergence checking

### 5. Analysis (`analysis.py`)
- Data analysis dashboard

### 6. Logging (`logger.py`)
- Standardized logging setup
- Color-coded output levels
- Consistent formatting

## Dependencies
- numpy: Numerical computations
- pandas: Data management
- scipy: Statistical operations
- dask: Delayed computations
- sklearn: Nearest neighbor calculations
- numpy_groupies: Group operations

