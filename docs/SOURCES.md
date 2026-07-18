# Source and validation boundary

The canonical JSON records source-level provenance in its top-level `sources` field. The current scientific checks use:

- the [Bilbao Crystallographic Server](https://www.cryst.ehu.es/) point-group retrieval tools for the displayed $D_{4h}$ and $D_{6h}$ operation order, coordinate action, and matrices;
- the user-supplied layer-group tables for LG1–LG80 point-group labels, $R^+/R^-$ membership, and the documented embeddings;
- exact matrix multiplication to generate and verify every $D_{4h}$ and $D_{6h}$ Cayley-table cell;
- the official [spglib Python API](https://spglib.readthedocs.io/en/stable/api/python-api.html) as the planned independent structure-symmetry integration boundary.

The repository does not redistribute the supplied reference PDF or screenshots. Those files are verification inputs, not package assets. The catalog does not claim that a matrix list alone reproduces Bilbao, International Tables, or spglib setting choices. Basis, origin, translations, tolerances, and settings remain part of any downstream comparison.

For Bilbao citation guidance and service scope, see the server's [About page](https://cryst.ehu.es/wiki/index.php/About_the_Bilbao_Crystallographic_Server). Access date for the public server links in this release candidate: 2026-07-17.
