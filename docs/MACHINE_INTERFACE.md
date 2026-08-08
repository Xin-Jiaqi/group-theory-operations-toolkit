# Machine interface and integration contract

## Stable public surface

```python
from group_theory_operations import (
    load_database,
    operation_record,
    multiply_operations,
    apply_fractional_operation,
    quadratic_field_representation,
    load_quadratic_field_catalog,
    load_point_group_registry,
    get_crystallographic_point_group,
    point_group_operations,
    response_tensor_basis,
    load_optical_response_catalog,
    load_magnetic_layer_group_registry,
    get_magnetic_layer_group,
    magnetic_layer_tensor_basis,
    magnetic_layer_response_tensor_basis,
)

database = load_database()
rotation = operation_record(database, "4+_001", "tetragonal_D4h")
product = multiply_operations(database, "tetragonal_D4h", "4+_001", "2_100")
transformed = apply_fractional_operation(structure, rotation)
fields = quadratic_field_representation(rotation.matrix_cartesian)
all_fields = load_quadratic_field_catalog()
registry = load_point_group_registry()
point_group = get_crystallographic_point_group("4mm")
operations = point_group_operations(point_group.number)
shift_basis = response_tensor_basis("4mm", "shift_current")
all_invariants = load_optical_response_catalog()
magnetic_layers = load_magnetic_layer_group_registry()
magnetic_layer = get_magnetic_layer_group("6.5.25", magnetic_layers)
odd_map = magnetic_layer_tensor_basis(
    magnetic_layer.og_number,
    "polar_vector",
    "symmetric_quadratic",
    input_time_parity="odd",
    registry=magnetic_layers,
)
shg_odd = magnetic_layer_response_tensor_basis(
    magnetic_layer.og_number, "shg_odd", registry=magnetic_layers
)
```

`OperationRecord` is frozen and exposes explicit fractional and Cartesian matrices. The original JSON remains the canonical data source, while `validate_database` provides a reusable schema/group-consistency gate. The CLI supports JSON output for downstream scripts:

```bash
group-ops list --json
group-ops show '4^+_{001}' --family tetragonal_D4h --json
group-ops multiply 4+_001 2_100 --family tetragonal_D4h --json
group-ops field-representation m_100 --family tetragonal_D4h --json
group-ops point-groups 4mm --json
group-ops invariants 4mm shift_current --json
group-ops magnetic-layer-groups 6.5.25 --json
group-ops magnetic-layer-responses 6.5.25 shg_odd --json
group-ops validate
```

## Quadratic-field representations

`quadratic_field_representation()` returns the two induced actions needed by
second-order optical-response workflows. Its default symmetric basis is

`(|Ex|², |Ey|², |Ez|², Ex Ey* + Ey Ex*, Ex Ez* + Ez Ex*, Ey Ez* + Ez Ey*)`.

The antisymmetric basis is the real axial vector `h = i E × E*`, and therefore
`matrix_antisymmetric = det(D) D`. The public functions derive both matrices
from the Cartesian orthogonal matrix; callers must not apply this axial formula
to the non-orthonormal fractional-coordinate representation.

`data/quadratic_field_representations.json` is a generated, versioned
all-operation view. It includes the source-catalog SHA-256 and follows
`schema/quadratic-field-representations-v1.schema.json`. Regenerate it with
`python scripts/generate_quadratic_field_representations.py`; the operation
catalog remains the sole primary source.

## Crystallographic point groups and optical invariants

`data/crystallographic_point_groups.json` registers the standard 32 point groups.
Each entry fixes its standard number, Hermann-Mauguin and Schoenflies symbols,
crystal system, host family, setting, generators, closure operations, and
centrosymmetric/polar/chiral flags. Operation names resolve through the primary
operation catalog; the registry does not duplicate matrices.

`equivariant_map_basis(A, B)` is the generic solver for maps satisfying
`A(R) T = T B(R)`. `response_tensor_basis(point_group, response)` applies it to
three spatial response contracts: `shift_current` and `shg` use the 3x6 map
`D <- M_+`; `circular_injection_current` uses the 3x3 map `D <- M_-`.
The returned deterministic basis spans every spatially allowed tensor. It does
not impose time reversal, frequency, resonance, units, or microscopic physics.

`data/optical_response_invariants.json` freezes the all-group results and binds
both registry and operation-catalog SHA-256 values. Regenerate the point-group
registry first, then the invariant catalog. Tests verify the closure and every
basis vector under every registered operation.

## Magnetic layer groups

`data/magnetic_layer_groups.json` registers all 528 magnetic layer groups by
global number, three-part OG number, type I–IV, parent layer group and
corresponding magnetic space group. Each record contains an explicit finite
magnetic point co-group `(R, time_reversal)` in both fractional and Cartesian
coordinates. Type-IV anti-translations are retained separately; ordinary
translations are intentionally omitted because this interface targets
homogeneous (`q=0`) tensor constraints, not affine action on atomic sites.

`magnetic_layer_tensor_basis()` handles arbitrary declared spatial spaces and
time parities. `magnetic_layer_response_tensor_basis()` provides the six named
shift-current, injection-current and SHG sectors. Exact affine magnetic-layer
operations require a future schema that also fixes origin and translation
conventions.

## Stacking-ferroelectricity boundary

`layer_group_polarization()` computes the common fixed space of the complete
Seitz-operation closure and classifies it as IP, OP, CP, or NP.
`stacking_rotation_cosets()` partitions a compatible Bravais-lattice point
group into left cosets of the monolayer point group; unlike a quotient-group
construction, the subgroup need not be normal. Both groups must use the same
coordinate embedding.

`equivalent_interface_translation()` implements the symmorphic relation
$R^{\pm}\tau_p=\pm\tau_q$. `preserves_recursive_stacking_step()` implements
the recursive multilayer condition: $R^+$ fixes both interface translations,
whereas $R^-$ exchanges them, modulo primitive translations and the optional
centered translation. These functions do not establish energetic degeneracy;
that premise must come from calculation or experiment before the returned
relations are interpreted as ferroelectric switching paths.

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

`schema/group-operations-v1.schema.json` describes the primary operation JSON. The quadratic-field, crystallographic-point-group, magnetic-layer-group and optical-invariant artifacts each have an independent v1 schema. Additive optional fields are backward compatible within each schema. Removing or renaming a field, changing matrix or multiplication semantics, or changing a stable operation name requires a new schema version and a migration note. JSON Schema checks shape and primitive types; scientific consistency remains covered by closure, spglib-signature, representation and equivariance tests.
