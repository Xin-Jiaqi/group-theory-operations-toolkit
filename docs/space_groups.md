# Crystallographic Space Groups (1–230)

The registry `data/crystallographic_space_groups.json` covers all 230
crystallographic space groups with full Seitz generators per Hall setting.
It closes the documented "non-zero translations" gap of the point-operation
catalog: space-group operations are genuine Seitz pairs $(R \mid \mathbf t)$
acting on fractional coordinates of the conventional cell.

## Contents

Every entry is keyed by ITA number and carries:

- international (short and full) and Schoenflies symbols;
- the parent point group (number, HM, Schoenflies) from the 32-point-group
  registry, with the crystal system inherited from it;
- the conventional-cell centering and the symmorphic flag;
- every Hall setting (all 530 Hall numbers are covered), each with its
  full operation count and a compact Seitz generating set.

The **primary** setting is the lowest Hall number whose operations
round-trip: applying them to three generic positions on a conventional
lattice of the right crystal system and asking spglib to identify the
resulting structure returns the same ITA number.

## Conventions

- Seitz action: $x' = R x + \mathbf t$ in fractional coordinates, with
  translations modulo the lattice.
- Generators are a deterministic greedy reduction; their closure
  reproduces the full operation count of the setting.
- The symmorphic flag uses the fixed-point criterion modulo the cell
  lattice: the group is symmorphic iff an origin shift makes every
  non-identity rotation pure.  Under this criterion exactly 73 of the
  230 groups are symmorphic, including $R3$ (146) and $R3m$ (160) in
  their rhombohedral-axis settings and excluding $Fd\bar 3m$ (227),
  whose diamond glide is essential.

## CLI

```bash
group-ops space-groups                 # table of all 230
group-ops space-groups 227             # one group
group-ops space-groups 227 --json      # machine-readable record
```

## Verification

The registry is machine-generated from the spglib database (BSD-3-Clause)
by `scripts/generate_crystallographic_space_groups.py`, which itself fails
unless the ITA crystal-system counts, the 73 symmorphic count, and the
32-point-group coverage hold.  `tests/test_space_groups.py` re-verifies the
published file from several independent directions:

1. structural invariants (counts, coverage, uniqueness, spot checks);
2. spglib cross-check (every stored label, Hall symbol, and operation
   count against `get_spacegroup_type` / `get_symmetry_from_database`);
3. generator closure (the stored generators close to the declared
   operation count for every one of the 530 settings);
4. ASE cross-check (HM symbols and operation-count divisibility against
   ASE's independent database; skipped when ASE is not installed);
5. gemmi cross-check (crystal system, operation count, and the symmorphic
   flag against gemmi's independent tables, which agree on all 230 groups;
   skipped when gemmi is not installed);
6. Wikipedia fixture (crystal system per ITA number transcribed from the
   Wikipedia list of space groups, agreeing on all 230);
7. end-to-end round trip (three generic positions per space group are
   identified by spglib as the same ITA number).

## Scope notes

- Operations are expressed in the conventional cell of each Hall setting;
  the primitive cell can be obtained by quotienting the identity-rotation
  translations (centering vectors) out of the operation set.
- The registry records the space groups themselves.  Wyckoff position
  multiplicities and site-symmetry tables are not included.
- The data is generated from spglib 2.x; regenerate after upgrading spglib
  and inspect the diff before committing.
