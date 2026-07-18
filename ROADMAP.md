# Release roadmap

## 0.1.0 candidate — packaged catalog

- [x] Keep one canonical JSON data source.
- [x] Add frozen operation records and JSON-capable CLI queries.
- [x] Add reusable catalog validation, including orthogonality and multiplication consistency.
- [x] Delegate structure I/O to `materials-structure-core`.
- [x] Verify operation composition against structure-coordinate transforms.
- [x] Reject incompatible PBC and Selective-dynamics transforms.
- [x] Build and inspect wheel/sdist artifacts.

## Stable-release gates

- [x] License repository-authored code, data, and documentation under BSD-3-Clause.
- [ ] Define a versioned JSON Schema and compatibility policy.
- [ ] Expand the registry toward all 32 crystallographic point groups with independent generator/order/closure tests.
- [ ] Cross-check structure classifications through a pinned spglib fixture suite.
- [ ] Define affine Seitz operations before accepting nonzero translations or shifted rotation centers.

The 32-point-group expansion is a data-curation project, not a name-only checklist. A group enters the registry only with sourced operations, a fixed basis/setting, closure tests, and a documented mapping to external conventions.
