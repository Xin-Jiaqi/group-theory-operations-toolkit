# Changelog

## 0.17.0 - 2026-09-27

- Add ordinary axial-spin/polar-current Edelstein tensor invariants for all
  crystallographic point groups.
- Add symmetry-only bilayer LEE classification with layer-preserving and
  layer-exchanging operation subsets, including component-wise Type-I,
  Type-II and Both results.
- Extend concrete-structure response screening to include Edelstein tensors.
- Correct stale README version and response-count descriptions.

## 0.16.0 - 2026-09-13

- Derive slab lattice point operations from the actual metric using materials-structure-core geometry.
- Classify polarization fixed spaces relative to a caller-supplied Cartesian layer normal.
- Preserve existing catalog APIs and base-only installation behavior.
- Add regression checks for all five 2D Bravais types, integer rebasing and rotations.
- Require materials-structure-core 0.0.3 in the optional structure extra.
