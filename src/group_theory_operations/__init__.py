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
from .representations import (
    QuadraticFieldRepresentation,
    antisymmetric_field_matrix,
    determinant3,
    load_quadratic_field_catalog,
    quadratic_field_representation,
    symmetric_field_matrix,
)

__all__ = [
    "GroupDataError",
    "OperationRecord",
    "QuadraticFieldRepresentation",
    "antisymmetric_field_matrix",
    "apply_fractional_operation",
    "canonical_name",
    "determinant3",
    "family_data",
    "find_operations",
    "get_layer_group",
    "get_operation",
    "get_point_group",
    "iter_operations",
    "load_database",
    "load_quadratic_field_catalog",
    "load_schema",
    "matrix_for",
    "multiply_operations",
    "operation_record",
    "quadratic_field_representation",
    "symmetric_field_matrix",
    "validate_database",
]

__version__ = "0.2.0"
