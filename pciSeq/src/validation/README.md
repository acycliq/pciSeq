# pciSeq Validation Module

This module handles validation and configuration management for the pciSeq pipeline.

## Files and Their Functions

### 1. `input_validation.py`
- Primary input validation module using `ValidatedInputs` dataclass
- Validates and processes core data structures:
  * Spots DataFrame (Gene, x, y coordinates)
  * COO matrix for segmentation masks
  * Single cell data (optional)

### 2. `config_manager.py`
- Manages configuration settings through `ConfigManager` dataclass
- Provides type hints for configuration parameters
