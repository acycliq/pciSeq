"""
Runtime Validation Module for pciSeq

This module handles runtime validations during algorithm execution,
including matrix dimensions, probability distributions, and convergence conditions.
"""

import numpy as np
import logging
from typing import Dict, Any, Optional, Tuple, Union
from dataclasses import dataclass

runtime_logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Container for validation results"""
    is_valid: bool
    message: str
    details: Optional[Dict[str, Any]] = None


def validate_probability_distribution(
        prob_matrix: np.ndarray,
        axis: int = -1,
        tolerance: float = 1e-10
) -> ValidationResult:
    """
    Validate probability distribution sums to 1 along specified axis.

    Args:
        prob_matrix: Array of probabilities
        axis: Axis along which probabilities should sum to 1
        tolerance: Maximum allowed deviation from 1

    Returns:
        ValidationResult with validation status and details
    """
    try:
        sums = np.sum(prob_matrix, axis=axis)
        is_valid = np.allclose(sums, 1.0, atol=tolerance)

        return ValidationResult(
            is_valid=is_valid,
            message="Probability distribution valid" if is_valid else "Invalid probability distribution",
            details={
                "max_deviation": float(np.max(np.abs(sums - 1.0))),
                "shape": prob_matrix.shape
            }
        )
    except Exception as e:
        runtime_logger.error(f"Probability validation failed: {str(e)}")
        raise


def validate_matrix_dimensions(
        *matrices: np.ndarray,
        expected_shapes: Dict[str, Tuple[int, ...]]
) -> ValidationResult:
    """
    Validate matrix dimensions match expected shapes.

    Args:
        *matrices: Arrays to validate
        expected_shapes: Dictionary mapping matrix names to expected shapes

    Returns:
        ValidationResult with validation status and details
    """
    try:
        actual_shapes = {name: mat.shape for name, mat in zip(expected_shapes.keys(), matrices)}
        mismatches = {
            name: (actual, expected)
            for name, (actual, expected) in zip(
                actual_shapes.keys(),
                zip(actual_shapes.values(), expected_shapes.values())
            )
            if actual != expected
        }

        is_valid = len(mismatches) == 0
        return ValidationResult(
            is_valid=is_valid,
            message="Matrix dimensions valid" if is_valid else "Invalid matrix dimensions",
            details={"mismatches": mismatches} if mismatches else None
        )
    except Exception as e:
        runtime_logger.error(f"Matrix dimension validation failed: {str(e)}")
        raise


def validate_convergence_state(
        delta: float,
        tolerance: float,
        iteration: int,
        max_iterations: int
) -> ValidationResult:
    """
    Validate convergence state of the algorithm.

    Args:
        delta: Change in values between iterations
        tolerance: Convergence tolerance threshold
        iteration: Current iteration number
        max_iterations: Maximum allowed iterations

    Returns:
        ValidationResult with validation status and details
    """
    try:
        is_valid = (delta < tolerance) or (iteration < max_iterations)

        return ValidationResult(
            is_valid=is_valid,
            message="Convergence state valid" if is_valid else "Invalid convergence state",
            details={
                "delta": delta,
                "tolerance": tolerance,
                "iteration": iteration,
                "max_iterations": max_iterations
            }
        )
    except Exception as e:
        runtime_logger.error(f"Convergence validation failed: {str(e)}")
        raise