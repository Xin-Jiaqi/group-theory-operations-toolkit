"""Symmetry kernels for bilayer and multilayer stacking ferroelectricity.

The routines in this module implement reusable group-theory statements rather
than material-specific energy models:

* allowed polarization is the common fixed space of all point operations;
* inequivalent bilayer orientations are left cosets of the monolayer point
  group in the Bravais-lattice point group;
* symmorphic, energetically equivalent interfaces obey
  ``R @ tau_1 = sign(R_zz) * tau_2``;
* a recursive multilayer step preserves an ``R+`` operation when both
  interface translations are fixed, and an ``R-`` operation when they are
  exchanged.

These statements follow the theory developed in Phys. Rev. Lett. 130, 146801
(2023), Phys. Rev. B 111, 224102 (2025), and Phys. Rev. B 113, 075310 (2026).
They determine symmetry allowance and switching relations.  They do not
predict interface energies, barriers, polarization magnitudes, or whether two
interfaces are energetically degenerate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Sequence, TypeVar

import numpy as np

from .catalog import GroupDataError, OperationRecord
from .layer_groups import layer_group_operations
from .point_groups import point_group_operations
from .seitz import SeitzOp


Matrix3 = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]
Vector2 = tuple[float, float]
Vector3 = tuple[float, float, float]

BRAVAIS_LATTICE_POINT_GROUPS = {
    "oP": "2/m",
    "rP": "mmm",
    "rC": "mmm",
    "sP": "4/mmm",
    "hP": "6/mmm",
}

_TOLERANCE = 1.0e-9
_OperationT = TypeVar("_OperationT")


@dataclass(frozen=True, slots=True)
class PolarizationSpace:
    """Common fixed space of a set of point operations.

    ``basis`` contains orthonormal column-space vectors expressed in the input
    coordinate system.  For non-Cartesian fractional coordinates the basis is
    algebraically valid, but its Euclidean normalization has no direct metric
    meaning.  ``polar_type`` is one of ``IP``, ``OP``, ``CP``, or ``NP``.
    """

    polar_type: str
    dimension: int
    in_plane_dimension: int
    out_of_plane_dimension: int
    basis: tuple[Vector3, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "polar_type": self.polar_type,
            "dimension": self.dimension,
            "in_plane_dimension": self.in_plane_dimension,
            "out_of_plane_dimension": self.out_of_plane_dimension,
            "basis": [list(vector) for vector in self.basis],
        }


@dataclass(frozen=True, slots=True)
class MatrixCoset:
    """One left coset and its selected representative."""

    representative: Matrix3
    members: tuple[Matrix3, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "representative": [list(row) for row in self.representative],
            "members": [
                [list(row) for row in matrix]
                for matrix in self.members
            ],
        }


@dataclass(frozen=True, slots=True)
class PolarizationSwitch:
    """Change of a polar vector under one Cartesian point operation."""

    initial: Vector3
    transformed: Vector3
    out_of_plane: str
    in_plane_angle_degrees: float | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "initial": list(self.initial),
            "transformed": list(self.transformed),
            "out_of_plane": self.out_of_plane,
            "in_plane_angle_degrees": self.in_plane_angle_degrees,
        }


def _matrix(operation: Any, *, coordinate: str = "fractional") -> np.ndarray:
    matrix: Any
    if isinstance(operation, SeitzOp):
        matrix = operation.rotation
    elif isinstance(operation, OperationRecord):
        matrix = operation.matrix(coordinate)
    else:
        matrix = operation
    result = np.asarray(matrix, dtype=np.float64)
    if result.shape != (3, 3) or not np.all(np.isfinite(result)):
        raise GroupDataError("point operation must be a finite 3x3 matrix")
    return result


def _matrix_tuple(matrix: np.ndarray) -> Matrix3:
    cleaned = np.asarray(matrix, dtype=np.float64).copy()
    cleaned[np.abs(cleaned) < _TOLERANCE] = 0.0
    rounded = np.rint(cleaned)
    cleaned[np.abs(cleaned - rounded) < _TOLERANCE] = rounded[
        np.abs(cleaned - rounded) < _TOLERANCE
    ]
    return tuple(tuple(float(value) for value in row) for row in cleaned)  # type: ignore[return-value]


def _matrix_key(matrix: np.ndarray) -> tuple[float, ...]:
    return tuple(float(value) for value in np.round(matrix, decimals=10).ravel())


def _validate_matrix_group(
    matrices: tuple[np.ndarray, ...], *, label: str
) -> dict[tuple[float, ...], np.ndarray]:
    by_key = {_matrix_key(matrix): matrix for matrix in matrices}
    if len(by_key) != len(matrices):
        raise GroupDataError(f"{label} operations must be unique")
    identity_key = _matrix_key(np.eye(3))
    if identity_key not in by_key:
        raise GroupDataError(f"{label} must contain the identity operation")
    for left in matrices:
        for right in matrices:
            if _matrix_key(left @ right) not in by_key:
                raise GroupDataError(f"supplied {label} operations are not closed")
    return by_key


def _in_plane_vector(translation: Any) -> np.ndarray:
    vector = np.asarray(translation, dtype=np.float64)
    if vector.shape == (3,):
        if abs(float(vector[2])) > _TOLERANCE:
            raise GroupDataError("interface translation must have zero z component")
        vector = vector[:2]
    if vector.shape != (2,) or not np.all(np.isfinite(vector)):
        raise GroupDataError("interface translation must be a finite pair or triplet")
    return vector


def _layer_matrix(operation: Any) -> tuple[np.ndarray, int]:
    matrix = _matrix(operation)
    if (
        np.max(np.abs(matrix[:2, 2])) > _TOLERANCE
        or np.max(np.abs(matrix[2, :2])) > _TOLERANCE
    ):
        raise GroupDataError("operation mixes the periodic plane and layer normal")
    sign = int(round(float(matrix[2, 2])))
    if sign not in {-1, 1} or abs(float(matrix[2, 2]) - sign) > _TOLERANCE:
        raise GroupDataError("operation must preserve or reverse the layer normal")
    return matrix, sign


def _canonical_translation(
    translation: Any, *, centering: str = "P"
) -> Vector2:
    vector = _in_plane_vector(translation) % 1.0
    vector[np.abs(vector) < _TOLERANCE] = 0.0
    if centering == "P":
        candidates = [vector]
    elif centering == "C":
        candidates = [vector, (vector - np.asarray((0.5, 0.5))) % 1.0]
    else:
        raise GroupDataError("centering must be 'P' or 'C'")
    canonical = min(
        candidates,
        key=lambda item: (round(float(item[0]), 10), round(float(item[1]), 10)),
    )
    canonical[np.abs(canonical) < _TOLERANCE] = 0.0
    return float(canonical[0]), float(canonical[1])


def _translations_equivalent(
    left: Any, right: Any, *, centering: str
) -> bool:
    return bool(
        np.allclose(
            _canonical_translation(left, centering=centering),
            _canonical_translation(right, centering=centering),
            atol=_TOLERANCE,
            rtol=0.0,
        )
    )


def polarization_space(operations: Iterable[Any]) -> PolarizationSpace:
    """Return the polar-vector subspace fixed by every supplied operation."""

    matrices = tuple(_matrix(operation) for operation in operations)
    if not matrices:
        raise GroupDataError("at least one point operation is required")
    constraints = np.vstack([matrix - np.eye(3) for matrix in matrices])
    _, singular_values, right = np.linalg.svd(constraints, full_matrices=True)
    rank = int(np.count_nonzero(singular_values > _TOLERANCE))
    raw_basis = right[rank:].T
    raw_basis[np.abs(raw_basis) < _TOLERANCE] = 0.0
    basis = tuple(
        tuple(float(value) for value in raw_basis[:, index])
        for index in range(raw_basis.shape[1])
    )

    in_plane_dimension = int(
        np.linalg.matrix_rank(raw_basis[:2, :], tol=_TOLERANCE)
    )
    out_of_plane_dimension = int(
        np.linalg.matrix_rank(raw_basis[2:3, :], tol=_TOLERANCE)
    )
    if raw_basis.shape[1] == 0:
        polar_type = "NP"
    elif out_of_plane_dimension == 0:
        polar_type = "IP"
    elif in_plane_dimension == 0:
        polar_type = "OP"
    else:
        polar_type = "CP"
    return PolarizationSpace(
        polar_type=polar_type,
        dimension=raw_basis.shape[1],
        in_plane_dimension=in_plane_dimension,
        out_of_plane_dimension=out_of_plane_dimension,
        basis=basis,  # type: ignore[arg-type]
    )


def layer_group_polarization(
    identifier: int | str, *, layer_hall_number: int | None = None
) -> PolarizationSpace:
    """Return the IP/OP/CP/NP classification of one crystallographic layer group."""

    return polarization_space(
        operation.rotation
        for operation in layer_group_operations(
            identifier, layer_hall_number=layer_hall_number
        )
    )


def partition_left_cosets(
    group_operations: Iterable[Any], subgroup_operations: Iterable[Any]
) -> tuple[MatrixCoset, ...]:
    """Partition a finite matrix group into left cosets ``g H``.

    Both sets must use the same coordinate embedding.  The function validates
    that every supplied subgroup operation is contained in the parent group and
    that every computed product remains in it.
    """

    group = tuple(_matrix(operation) for operation in group_operations)
    subgroup = tuple(_matrix(operation) for operation in subgroup_operations)
    if not group or not subgroup:
        raise GroupDataError("group and subgroup must both be non-empty")
    group_by_key = _validate_matrix_group(group, label="parent group")
    subgroup_by_key = _validate_matrix_group(subgroup, label="subgroup")
    subgroup_keys = set(subgroup_by_key)
    if not subgroup_keys.issubset(group_by_key):
        raise GroupDataError(
            "subgroup is not embedded in the parent group in this coordinate basis"
        )

    covered: set[tuple[float, ...]] = set()
    cosets: list[MatrixCoset] = []
    for representative in group:
        representative_key = _matrix_key(representative)
        if representative_key in covered:
            continue
        member_keys = tuple(_matrix_key(representative @ item) for item in subgroup)
        if len(set(member_keys)) != len(subgroup):
            raise GroupDataError("left-coset multiplication produced duplicates")
        missing = [key for key in member_keys if key not in group_by_key]
        if missing:  # pragma: no cover - guarded by group/subgroup validation
            raise GroupDataError("left-coset product is outside the parent group")
        covered.update(member_keys)
        cosets.append(
            MatrixCoset(
                representative=_matrix_tuple(representative),
                members=tuple(
                    _matrix_tuple(group_by_key[key]) for key in member_keys
                ),
            )
        )
    if len(covered) != len(group):
        raise GroupDataError("left cosets do not partition the supplied parent group")
    return tuple(cosets)


def bravais_lattice_operations(lattice_type: str) -> tuple[OperationRecord, ...]:
    """Return point operations for one of the five two-dimensional Bravais lattices."""

    try:
        point_group = BRAVAIS_LATTICE_POINT_GROUPS[lattice_type]
    except KeyError as exc:
        choices = ", ".join(BRAVAIS_LATTICE_POINT_GROUPS)
        raise GroupDataError(
            f"unknown 2D Bravais lattice {lattice_type!r}; choices: {choices}"
        ) from exc
    return point_group_operations(point_group)


def stacking_rotation_cosets(
    monolayer_point_group: str | int, lattice_type: str
) -> tuple[MatrixCoset, ...]:
    """Return inequivalent same-cell bilayer rotation classes.

    The standard point-group embedding must be a subgroup of the requested
    Bravais-lattice point group.  If not, callers must supply the physically
    correct embedding directly to :func:`partition_left_cosets`.
    """

    lattice = bravais_lattice_operations(lattice_type)
    monolayer = point_group_operations(monolayer_point_group)
    return partition_left_cosets(
        (operation.matrix_fractional for operation in lattice),
        (operation.matrix_fractional for operation in monolayer),
    )


def equivalent_interface_translation(
    translation: Sequence[float],
    operation: Any,
    *,
    centering: str = "P",
) -> Vector2:
    """Map one interface translation to a symmetry-equivalent translation.

    For a layer-preserving operation ``R+`` this returns ``R @ tau``.  For a
    layer-exchanging operation ``R-`` it returns ``-R @ tau``, as required by
    the bilayer relation ``R^± tau_p = ± tau_q``.
    """

    matrix, sign = _layer_matrix(operation)
    vector = sign * (matrix[:2, :2] @ _in_plane_vector(translation))
    return _canonical_translation(vector, centering=centering)


def equivalent_interface_orbit(
    translation: Sequence[float],
    operations: Iterable[Any],
    *,
    centering: str = "P",
) -> tuple[Vector2, ...]:
    """Return the distinct symmetry-equivalent translations of one interface."""

    orbit: dict[tuple[float, float], Vector2] = {}
    for operation in operations:
        vector = equivalent_interface_translation(
            translation, operation, centering=centering
        )
        key = (round(float(vector[0]), 10), round(float(vector[1]), 10))
        orbit[key] = vector
    if not orbit:
        raise GroupDataError("at least one point operation is required")
    return tuple(orbit[key] for key in sorted(orbit))


def preserves_recursive_stacking_step(
    operation: Any,
    translation_1: Sequence[float],
    translation_2: Sequence[float],
    *,
    centering: str = "P",
) -> bool:
    """Test the recursive multilayer preservation criterion.

    ``R+`` must fix both interfacial translations.  ``R-`` must exchange them.
    Equality is modulo primitive lattice translations and, for ``centering=C``,
    modulo the centered translation ``(1/2, 1/2)``.
    """

    matrix, sign = _layer_matrix(operation)
    first = _in_plane_vector(translation_1)
    second = _in_plane_vector(translation_2)
    mapped_first = matrix[:2, :2] @ first
    mapped_second = matrix[:2, :2] @ second
    if sign == 1:
        return _translations_equivalent(
            mapped_first, first, centering=centering
        ) and _translations_equivalent(
            mapped_second, second, centering=centering
        )
    return _translations_equivalent(
        mapped_first, second, centering=centering
    ) and _translations_equivalent(
        mapped_second, first, centering=centering
    )


def preserved_recursive_stacking_operations(
    operations: Iterable[_OperationT],
    translation_1: Sequence[float],
    translation_2: Sequence[float],
    *,
    centering: str = "P",
) -> tuple[_OperationT, ...]:
    """Filter point operations by the recursive multilayer criterion."""

    return tuple(
        operation
        for operation in operations
        if preserves_recursive_stacking_step(
            operation,
            translation_1,
            translation_2,
            centering=centering,
        )
    )


def polarization_switch(
    operation_cartesian: Any, polarization: Sequence[float]
) -> PolarizationSwitch:
    """Describe OP reversal and the IP rotation angle under a Cartesian operation."""

    matrix = _matrix(operation_cartesian, coordinate="cartesian")
    vector = np.asarray(polarization, dtype=np.float64)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise GroupDataError("polarization must be a finite length-3 vector")
    if float(np.linalg.norm(vector)) < _TOLERANCE:
        raise GroupDataError("polarization must be non-zero")
    transformed = matrix @ vector
    if not math.isclose(
        float(np.linalg.norm(transformed)),
        float(np.linalg.norm(vector)),
        rel_tol=1.0e-8,
        abs_tol=1.0e-8,
    ):
        raise GroupDataError("operation_cartesian must preserve vector length")

    if abs(float(vector[2])) < _TOLERANCE:
        out_of_plane = "absent"
    elif float(vector[2] * transformed[2]) > 0:
        out_of_plane = "unchanged"
    elif float(vector[2] * transformed[2]) < 0:
        out_of_plane = "reversed"
    else:
        out_of_plane = "changed"

    in_plane_norm = float(np.linalg.norm(vector[:2]))
    transformed_in_plane_norm = float(np.linalg.norm(transformed[:2]))
    if in_plane_norm < _TOLERANCE or transformed_in_plane_norm < _TOLERANCE:
        angle = None
    else:
        cosine = float(
            np.dot(vector[:2], transformed[:2])
            / (in_plane_norm * transformed_in_plane_norm)
        )
        angle = float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        if abs(angle) < 1.0e-8:
            angle = 0.0
        elif abs(angle - round(angle)) < 1.0e-8:
            angle = float(round(angle))
    return PolarizationSwitch(
        initial=tuple(float(value) for value in vector),  # type: ignore[arg-type]
        transformed=tuple(float(value) for value in transformed),  # type: ignore[arg-type]
        out_of_plane=out_of_plane,
        in_plane_angle_degrees=angle,
    )
