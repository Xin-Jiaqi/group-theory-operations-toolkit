"""Machine-readable crystallographic point-operation catalog."""

from .catalog import (
    GroupDataError,
    OperationRecord,
    canonical_name,
    family_data,
    find_operations,
    get_layer_group,
    get_operation,
    get_point_group,
    iter_operations,
    load_database,
    load_schema,
    matrix_for,
    multiply_operations,
    operation_record,
    validate_database,
)
from .structure import apply_fractional_operation

__all__ = [
    "GroupDataError",
    "OperationRecord",
    "apply_fractional_operation",
    "canonical_name",
    "family_data",
    "find_operations",
    "get_layer_group",
    "get_operation",
    "get_point_group",
    "iter_operations",
    "load_database",
    "load_schema",
    "matrix_for",
    "multiply_operations",
    "operation_record",
    "validate_database",
]

__version__ = "0.1.0"
