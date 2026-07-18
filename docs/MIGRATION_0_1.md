# Migration to the 0.1 candidate

## Preserved query behavior

The source-tree `group_tools.py` launcher still supports `list`, `show`, `group`, and `multiply`. Existing query functions such as `canonical_name`, `load_database`, `get_operation`, `matrix_for`, and `multiply_operations` remain importable from that compatibility module.

New integrations should import from `group_theory_operations` and use frozen `OperationRecord` values. `load_database()` now validates the catalog's core schema, referenced metadata, matrices, subgroup structure and multiplication data by default; callers intentionally inspecting a damaged catalog may pass `validate=False`. Descriptive top-level prose and citation records are not a general-purpose JSON Schema contract.

## Intentional breaking change

The handwritten `parse_poscar`, `write_poscar`, `transform_coordinates`, and `apply-poscar` path has been removed. It supported only a narrow VASP 5 subset and could silently mishandle scale/species/format semantics. Replace it with:

```bash
group-ops apply-structure input.vasp 4+_001 \
  --family tetragonal_D4h --output output.vasp
```

Install the sibling `materials-structure-core[io]>=0.0.2` candidate for this command. Programmatic consumers should call its `read_structure`/`write_structure` functions and this package's `apply_fractional_operation`.

## Semantics now enforced

- matrices act on column vectors;
- multiplication means row/left times column/right, so the right operation acts first;
- structure transforms are about the origin and keep the lattice fixed;
- incompatible PBC axes and unrepresentable Selective-dynamics mappings raise `GroupDataError`.
