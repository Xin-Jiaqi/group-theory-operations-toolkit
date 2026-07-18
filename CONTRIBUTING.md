# Contributing

- Treat `data/group_operations.json` as the only canonical data source.
- State coordinate basis, vector convention, setting/embedding, source page, access date, and tolerance for every scientific change.
- Add a failing regression test before changing an operation name, index, matrix, group membership, or multiplication result.
- Do not add nonzero translations to the point-operation schema; design and version an affine Seitz schema first.
- Do not parse POSCAR/CIF here. Extend the shared `materials-structure-core` adapter and consume its public contract.
- Do not commit external tables or structures without provenance and redistribution review.
- Require independent scientific review for catalog changes and a separate software review for API/schema changes.
