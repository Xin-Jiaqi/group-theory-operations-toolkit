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

## 0.4.0 — crystallographic space and layer groups

- [x] Register all 230 space groups with international/Schoenflies labels, parent point groups, crystal systems, centering and symmorphic flags.
- [x] Cover every one of the 530 Hall settings with Seitz generators and operation counts.
- [x] Round-trip-verify the primary setting of every space group through spglib identification.
- [x] Cross-check labels, Hall symbols, and operation counts against spglib, and HM symbols against ASE.
- [x] Fix the symmorphic flag by the lattice-aware fixed-point criterion (73 groups).
- [x] Register all 80 layer groups and all 116 layer Hall settings with Seitz generators.
- [x] Cross-link every layer group to the existing axial point-group catalog and verify every generated operation against spglib.
- [x] Fix Seitz inversion for non-orthogonal fractional-coordinate bases and declare NumPy as a core runtime dependency.

## 0.4.1 — closure and layer cross-check maintenance

- [x] Accept reusable collections and one-shot iterators in Seitz closure.
- [x] Cross-check the standard rotation set of every LG1–LG80 entry against the legacy point-operation catalog.

## 0.5.0 — magnetic point groups and time-parity tensors

- [x] Register all 122 magnetic point groups with standard three-part numbers and traditional Hermann–Mauguin symbols.
- [x] Generate type-I, gray, and black-white operation sets from the validated 32-point-group embeddings.
- [x] Verify all colored groups algebraically and cross-check their complete classification against spglib's 1651 magnetic space groups.
- [x] Add Python and CLI access to magnetic operations and real tensor maps with explicit input/output time parity.

## 0.6.0 — magnetic nonlinear-optical sectors

- [x] Separate normal/magnetic shift and injection current sectors by polarization and magnetic-domain character.
- [x] Add i-type/c-type SHG bases for all 122 magnetic point groups.
- [x] Add a complex antiunitary map solver with explicit conjugation and caller-defined right actions for frequency-channel permutations.
- [x] Regression-test gray, $PT$, and published $\bar3'm'$ selection rules.
- [ ] Add response-specific causal identities between positive- and negative-frequency susceptibilities only after the damping and Fourier conventions are selected by the consuming model.

## 0.7.0 — all 528 magnetic layer groups

- [x] Register all 528 magnetic layer groups with OG numbering and type I–IV classification.
- [x] Store finite magnetic point-co-group operations and preserve every type-IV anti-translation separately.
- [x] Link every record to its parent layer group and corresponding magnetic space group.
- [x] Cross-check numbering and correspondence against two independent published tables and spglib.
- [x] Solve the six magnetic nonlinear-optical sectors for every magnetic layer group.
- [x] Add versioned JSON/Schema artifacts, deterministic extraction/generation scripts, typed Python APIs and CLI queries.

## 0.8.0 — bilayer and multilayer stacking symmetry

- [x] Derive IP/OP/CP/NP classifications from the fixed polar-vector subspace and reproduce all 80 published layer-group classifications.
- [x] Replace factor-group assumptions by a left-coset partition that also accepts nonnormal monolayer subgroups.
- [x] Implement the symmorphic equivalent-interface relation for layer-preserving and layer-exchanging operations.
- [x] Implement the recursive multilayer preservation criterion for primitive and centered in-plane lattices.
- [x] Reproduce the BN AB/BA relation and the graphene ABC/$D_{3d}$ and ABA/$D_{3h}$ benchmarks.
- [x] Keep energetic degeneracy, switching barriers and polarization magnitudes outside the symmetry-only contract.

## Stable-release gates

- [x] License repository-authored code, data, and documentation under BSD-3-Clause.
- [x] Define a versioned JSON Schema and compatibility policy.
- [x] Expand the registry to all 32 crystallographic point groups with independent generator/order/closure tests.
- [x] Add invariant-tensor solvers on top of the stable $M_+$/$M_-$ interfaces.
- [x] Cross-check real-structure classifications in all seven crystal systems through fixtures pinned to spglib v2.5.0 source files and hashes.
- [ ] Define affine Seitz operations on concrete structures (origin choice, Wyckoff placement) before accepting them for structure transforms.

The 32-point-group expansion is a data-curation project, not a name-only checklist. A group enters the registry only with sourced operations, a fixed basis/setting, closure tests, and a documented mapping to external conventions.
