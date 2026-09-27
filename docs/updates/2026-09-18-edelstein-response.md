# Edelstein response update

## What changed

The toolkit supports the ordinary Edelstein response tensor through
`response_tensor_basis(point_group, "edelstein")` and the CLI command
`group-ops invariants <point-group> edelstein --json`.  It now also supports
the symmetry-only bilayer layer Edelstein effect (LEE) through
`bilayer_edelstein_point_groups(monolayer, bilayer)` and
`group-ops bilayer-edelstein <monolayer> <bilayer> --json`.

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
- The bilayer classifier partitions $R_B^+$ and $R_B^-$ operations, reports
  Type-I / Type-II / Both component sets, and reproduces the MoSSe Type-I,
  MoTe$_2$ Type-II, and $D_3\rightarrow C_i$ Both routes.

## Scope

The LEE result is a symmetry-allowance classification.  It does not calculate
response magnitudes, interlayer tunnelling, Fermi-level transport weights, or
the electrostatic layer detuning caused by a finite vertical field.
