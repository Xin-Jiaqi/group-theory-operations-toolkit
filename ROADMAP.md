# Release roadmap

## 0.1.0 candidate — packaged catalog

- [x] Keep one canonical JSON data source.
- [x] Add frozen operation records and JSON-capable CLI queries.
- [x] Add reusable catalog validation, including orthogonality and multiplication consistency.
- [x] Delegate structure I/O to `materials-structure-core`.
- [x] Verify operation composition against structure-coordinate transforms.
- [x] Reject incompatible PBC and Selective-dynamics transforms.
- [x] Build and inspect wheel/sdist artifacts.

## 0.2.0 candidate — quadratic-field representations

- [x] Derive $M_+(R)=\operatorname{Sym}^2D(R)$ in the explicit complex-field basis.
- [x] Derive $M_-(R)=\det[D(R)]D(R)$ for the axial-vector input space.
- [x] Publish generated 6×6 and 3×3 matrices for all 88 family records.
- [x] Add a versioned derived-data schema, Python API and JSON CLI query.
- [x] Verify direct complex-field action, cross products, orthogonality and representation homomorphisms.

## 0.3.0 candidate — crystallographic registry and optical invariants

- [x] Register all 32 crystallographic point groups with fixed settings, generators and closures.
- [x] Cross-check group order and element signatures against spglib's Hall database.
- [x] Solve spatial invariant bases for shift current, SHG and circular injection current.
- [x] Publish versioned JSON/Schema artifacts, Python APIs, CLI queries and readable dimension tables.
- [x] Verify every generated basis vector against every operation of every point group.

## Stable-release gates

- [x] License repository-authored code, data, and documentation under BSD-3-Clause.
- [x] Define a versioned JSON Schema and compatibility policy.
- [x] Expand the registry to all 32 crystallographic point groups with independent generator/order/closure tests.
- [x] Add invariant-tensor solvers on top of the stable $M_+$/$M_-$ interfaces.
- [ ] Cross-check structure classifications through a pinned spglib fixture suite.
- [ ] Define affine Seitz operations before accepting nonzero translations or shifted rotation centers.

The 32-point-group expansion is a data-curation project, not a name-only checklist. A group enters the registry only with sourced operations, a fixed basis/setting, closure tests, and a documented mapping to external conventions.
