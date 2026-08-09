"""Symmetry-allowed nonlinear responses for concrete crystal structures."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

from .invariants import ResponseSymmetryResult, screen_response_symmetry
from .structure_symmetry import StructureSymmetryContext, classify_structure_symmetry


@dataclass(frozen=True, slots=True)
class StructureResponseAnalysis:
    """Space-group classification and response selection for one structure.

    The response dimensions are obtained from the crystallographic point group
    associated with the detected space group.  They are independent of a rigid
    rotation of the crystal axes.  A nonzero dimension means only that the
    response is not forbidden by spatial symmetry; it does not determine a
    material-specific magnitude, spectrum, or experimental observability.
    """

    symmetry: StructureSymmetryContext
    responses: tuple[ResponseSymmetryResult, ...]

    @property
    def allowed_responses(self) -> tuple[ResponseSymmetryResult, ...]:
        """Return only response sectors with a nonzero invariant dimension."""

        return tuple(result for result in self.responses if result.allowed)

    def to_dict(self) -> dict[str, Any]:
        """Return a compact JSON-compatible physical summary."""

        space_group = self.symmetry.space_group
        hall_setting = self.symmetry.hall_setting
        return {
            "space_group": {
                "ita_number": space_group.ita_number,
                "international_short": space_group.international_short,
                "international_full": space_group.international_full,
                "schoenflies_symbol": space_group.schoenflies,
                "crystal_system": space_group.crystal_system,
            },
            "hall_setting": {
                "hall_number": hall_setting.hall_number,
                "hall_symbol": hall_setting.hall_symbol,
                "choice": hall_setting.choice,
                "centering": hall_setting.centering,
            },
            "point_group": {
                "number": space_group.point_group_number,
                "hm_symbol": space_group.point_group_hm,
                "schoenflies_symbol": space_group.point_group_schoenflies,
            },
            "input_structure": {
                "site_count": len(self.symmetry.equivalent_atoms),
                "equivalent_atom_orbit_count": len(
                    set(self.symmetry.equivalent_atoms)
                ),
                "operation_count": len(self.symmetry.input_operations),
                "wyckoff_letters": list(self.symmetry.wyckoff_letters),
                "site_symmetry_symbols": list(
                    self.symmetry.site_symmetry_symbols
                ),
            },
            "standard_setting": {
                "operation_count": len(self.symmetry.standard_operations),
                "transformation_matrix": [
                    list(row) for row in self.symmetry.transformation_matrix
                ],
                "origin_shift": list(self.symmetry.origin_shift),
            },
            "classification": {
                "symprec": self.symmetry.symprec,
                "angle_tolerance": self.symmetry.angle_tolerance,
                "backend": self.symmetry.backend_name,
                "backend_version": self.symmetry.backend_version,
            },
            "responses": [result.to_dict() for result in self.responses],
        }


def analyze_structure_responses(
    structure: Any,
    *,
    responses: Iterable[str] | None = None,
    allowed_only: bool = False,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
    database: Mapping[str, Any] | None = None,
    point_group_registry: Mapping[str, Any] | None = None,
    tolerance: float = 1.0e-8,
) -> StructureResponseAnalysis:
    """Classify a nonmagnetic crystal and screen its second-order responses.

    The supplied structure must be periodic along all three lattice directions.
    The spatial symmetry search is delegated to spglib and independently checked
    against the package Hall registry by :func:`classify_structure_symmetry`.
    The detected crystallographic point-group number is then used to evaluate
    shift current, electric-dipole SHG, and circular injection current through
    character inner products.

    Magnetic response sectors are intentionally excluded because atomic species
    and positions alone do not specify magnetic moments or a magnetic group.
    """

    symmetry = classify_structure_symmetry(
        structure,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    screened = screen_response_symmetry(
        "point_group",
        groups=(symmetry.space_group.point_group_number,),
        responses=responses,
        allowed_only=allowed_only,
        database=database,
        registry=point_group_registry,
        tolerance=tolerance,
    )
    return StructureResponseAnalysis(symmetry=symmetry, responses=screened)
