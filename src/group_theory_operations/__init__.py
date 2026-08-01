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
from .point_groups import (
    CrystallographicPointGroup,
    get_crystallographic_point_group,
    iter_crystallographic_point_groups,
    load_point_group_registry,
    point_group_operations,
)
from .invariants import (
    InvariantTensorBasis,
    canonical_response_name,
    equivariant_map_basis,
    load_optical_response_catalog,
    response_tensor_basis,
)

__all__ = [
    "GroupDataError",
    "CrystallographicPointGroup",
    "InvariantTensorBasis",
    "OperationRecord",
    "QuadraticFieldRepresentation",
    "antisymmetric_field_matrix",
    "apply_fractional_operation",
    "canonical_name",
    "canonical_response_name",
    "determinant3",
    "family_data",
    "find_operations",
    "get_crystallographic_point_group",
    "get_layer_group",
    "get_operation",
    "get_point_group",
    "iter_operations",
    "iter_crystallographic_point_groups",
    "load_database",
    "load_optical_response_catalog",
    "load_point_group_registry",
    "load_quadratic_field_catalog",
    "load_schema",
    "matrix_for",
    "multiply_operations",
    "operation_record",
    "point_group_operations",
    "quadratic_field_representation",
    "response_tensor_basis",
    "symmetric_field_matrix",
    "equivariant_map_basis",
    "validate_database",
]

__version__ = "0.3.0"
