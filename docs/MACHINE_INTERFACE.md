# Machine interface and integration contract

## Stable surface for the 0.1 candidate

```python
from group_theory_operations import (
    load_database,
    operation_record,
    multiply_operations,
    apply_fractional_operation,
)

database = load_database()
rotation = operation_record(database, "4+_001", "tetragonal_D4h")
product = multiply_operations(database, "tetragonal_D4h", "4+_001", "2_100")
transformed = apply_fractional_operation(structure, rotation)
```

`OperationRecord` is frozen and exposes explicit fractional and Cartesian matrices. The original JSON remains the canonical data source, while `validate_database` provides a reusable schema/group-consistency gate. The CLI supports JSON output for downstream scripts:

```bash
group-ops list --json
group-ops show '4^+_{001}' --family tetragonal_D4h --json
group-ops multiply 4+_001 2_100 --family tetragonal_D4h --json
group-ops validate
```

## Structure-core boundary

`apply_fractional_operation` implements the column-vector rule $f' = Df$ about the origin and keeps the lattice fixed. It accepts the public `materials-structure-core` `StructureRecord` contract. The bridge preserves site order and supported site flags.

The function rejects:

- operations whose fractional matrix does not preserve the supplied row-lattice metric, `D.T @ (L @ L.T) @ D = L @ L.T`;
- operations that mix periodic and non-periodic directions;
- operations that change the declared PBC-axis semantics;
- non-signed-permutation matrices when component-wise Selective-dynamics flags are present.

These are contract failures, not numerical warnings. A Boolean flag attached to one lattice-coordinate direction cannot generally be mapped exactly by a hexagonal fractional matrix that mixes two directions.

Family basis, setting and origin are part of the operation contract. For
example, a tetragonal fourfold matrix is rejected on an orthorhombic cell even
though it is still an invertible fractional-coordinate transform. This API is
therefore not an unrestricted algebraic coordinate transformer.

`group-ops apply-structure` delegates POSCAR/CIF parsing and writing to `materials-structure-core[io]>=0.0.2`. This repository no longer maintains its own POSCAR parser.

## Scientific convention

Matrices act on column vectors. Therefore `multiply_operations(left, right)` represents $D(left)D(right)$ and the right operation acts first. Integration tests verify that the table result and sequential structure transforms agree.

## Independent structure validation

[spglib's symmetry dataset API](https://spglib.readthedocs.io/en/stable/api/python-api.html#spglib.spglib.SpglibDataset)
is the intended independent structure-symmetry consumer/checker. This repository
does not claim to replace symmetry inference or standardization.
# Schema compatibility

`schema/group-operations-v1.schema.json` describes the public JSON shape for `schema_version: 1`. Additive optional fields are backward compatible within schema v1. Removing or renaming a field, changing matrix or multiplication semantics, or changing a stable operation name requires a new schema version and a migration note. The JSON Schema checks shape and primitive types; `validate_database()` remains authoritative for group closure, basis conversion, inverse, and multiplication consistency.
