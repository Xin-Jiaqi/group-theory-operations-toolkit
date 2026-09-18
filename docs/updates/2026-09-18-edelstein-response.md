# Edelstein response update

## What changed

The toolkit now supports the ordinary Edelstein response tensor through
`response_tensor_basis(point_group, "edelstein")` and the CLI command
`group-ops invariants <point-group> edelstein --json`.

It solves $S_i=\chi_{ij}j_j$, treating spin as an axial vector and current
as a polar vector.  The symmetry constraint is

$$
\chi=\det(R)R\chi R^T.
$$

## Checks and examples

- Inversion symmetry forces $\chi=0$.
- For $C_{4v}$ (`4mm`), the only basis is the Rashba form
  $\chi_{xy}=-\chi_{yx}$.
- Regression tests cover centrosymmetric, mirror, dihedral, and hexagonal
  representatives.

## Scope

This update covers ordinary homogeneous Edelstein selection rules.  The
layer-exchange and field-induced classification needed for bilayer LEE
(Type-I / Type-II) remains a separate future extension.
