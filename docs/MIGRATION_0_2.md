# Migration to the 0.2 candidate

Version 0.2 adds induced representations without changing operation names,
matrices, multiplication semantics, or the schema-v1 operation catalog.

New public interfaces:

- `symmetric_field_matrix(D, normalization="sum")`;
- `antisymmetric_field_matrix(D)`;
- `quadratic_field_representation(operation_or_matrix)`;
- `load_quadratic_field_catalog()`;
- `group-ops field-representation ...`.

The default six-dimensional basis uses the full off-diagonal sums. Historical
matrices in `docs/group_theory.md` used half of each sum. They are related by
`M_sum = S M_half S^-1`, where `S = diag(1, 1, 1, 2, 2, 2)`; no physics changes,
but code must not mix the two coordinate conventions.

The new `data/quadratic_field_representations.json` is derived from
`data/group_operations.json`. It is safe for machine consumption but is not a
second manually editable source.
