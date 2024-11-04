# pciSeq Validation Module

This module handles validation and configuration management for the pciSeq pipeline.

## Files and Their Functions

### 1. `input_validation.py`
- Primary input validation module using `ValidatedInputs` dataclass
- Validates and processes core data structures:
  * Spots DataFrame (Gene, x, y coordinates)
  * COO matrix for segmentation masks
  * Single cell data (optional)
- Handles data type casting and gene filtering

### 2. `config_manager.py`
- Manages configuration settings through `ConfigManager` dataclass
- Handles parameter validation and type checking:
  * InsideCellBonus settings
  * MisreadDensity configuration
  * cell_type_prior options
- Provides configuration conversion between dict and class formats

### 3. `runtime_validation.py`
- Handles runtime validations during algorithm execution
- Provides validation for:
  * Probability distributions (`validate_probability_distribution`)
  * Matrix dimensions (`validate_matrix_dimensions`)
  * Convergence conditions (`validate_convergence_state`)
