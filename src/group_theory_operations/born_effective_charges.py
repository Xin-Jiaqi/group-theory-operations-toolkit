"""Static Born-charge selection rules, atom mappings, and bilayer mode charges.

All tensors use the raw VASP OUTCAR convention: Z[k, i, j] has the
polarization/electric-field index i first and displacement/force index j
second.  Consequently F[k] = e * Z[k].T @ E, not e * Z[k] @ E.

The solver accepts generators or a full group, in orthonormal Cartesian
coordinates.  It does not assume Z is symmetric.  A site-group answer is
not a whole-crystal answer: use ``structure_born_space`` to add atom
permutations and the neutral-cell acoustic sum rule.

Only NumPy is required.  Named point groups use the toolkit's existing
registry lazily; no second runtime point-group database is introduced.
See docs/born-effective-charges.md for assumptions and derivations.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]
DEFAULT_TOL = 1e-10
CONVENTION = "Z[k,i,j]: i=polarization/field, j=displacement/force (raw VASP OUTCAR)"


def _tolerance(tol: float) -> float:
    if not np.isfinite(tol) or tol <= 0 or tol >= 1e-2:
        raise ValueError("tol must be finite, positive, and smaller than 1e-2")
    return float(tol)


def _real_array(value: ArrayLike, name: str) -> FloatArray:
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real (static adiabatic response)")
    try:
        array = np.asarray(value, dtype=float)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a real numerical array") from exc
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains a non-finite value")
    return array


def _rotations(operations: Iterable[Any], tol: float) -> list[FloatArray]:
    result: list[FloatArray] = []
    for index, operation in enumerate(operations):
        # Accept an existing toolkit OperationRecord as well as a raw matrix.
        value = getattr(operation, "matrix_cartesian", operation)
        matrix = _real_array(value, f"operation {index}")
        if matrix.shape != (3, 3):
            raise ValueError(f"operation {index} must have shape (3, 3)")
        if not np.allclose(matrix.T @ matrix, np.eye(3), atol=tol, rtol=0):
            raise ValueError("operations must be orthonormal Cartesian matrices; "
                             "convert fractional operations with R=L W L^-1 first")
        result.append(matrix)
    if not result:
        raise ValueError("operations cannot be empty; supply identity for C1")
    return result


@dataclass(frozen=True)
class BornChargeSpace:
    """An orthonormal basis of an allowed real tensor/vector subspace.

    ``basis.shape == (dimension, *tensor_shape)``.  Rows of the flattened
    basis are orthonormal in the Frobenius inner product.  Numerical basis
    vectors may mix independent parameters; compare projectors rather
    than SVD basis entries.  ``allowed_mask`` means *not symmetry-forced
    to zero*, not a prediction of a finite physical coefficient.
    """

    tensor_shape: tuple[int, ...]
    basis: FloatArray
    tolerance: float = DEFAULT_TOL

    def __post_init__(self) -> None:
        _tolerance(self.tolerance)
        if not self.tensor_shape or any(type(n) is not int or n <= 0 for n in self.tensor_shape):
            raise ValueError("tensor_shape must contain positive integers")
        array = _real_array(self.basis, "basis").copy()
        if array.ndim != len(self.tensor_shape) + 1 or array.shape[1:] != self.tensor_shape:
            raise ValueError("basis shape does not match tensor_shape")
        flat = array.reshape((array.shape[0], int(np.prod(self.tensor_shape))))
        if not np.allclose(flat @ flat.T, np.eye(len(flat)), atol=20*self.tolerance, rtol=0):
            raise ValueError("basis must be orthonormal")
        array.setflags(write=False)
        object.__setattr__(self, "basis", array)

    @property
    def dimension(self) -> int:
        return self.basis.shape[0]

    @property
    def allowed_mask(self) -> NDArray[np.bool_]:
        return np.linalg.norm(self.basis, axis=0) > self.tolerance

    @property
    def projector(self) -> FloatArray:
        flat = self.basis.reshape((self.dimension, int(np.prod(self.tensor_shape))))
        return flat.T @ flat

    def project(self, tensor: ArrayLike) -> FloatArray:
        """Orthogonally project a numerical tensor; no silent ASR correction."""
        array = _real_array(tensor, "tensor")
        if array.shape != self.tensor_shape:
            raise ValueError(f"tensor must have shape {self.tensor_shape}")
        return (self.projector @ array.ravel()).reshape(self.tensor_shape)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "convention": CONVENTION,
            "tensor_shape": list(self.tensor_shape),
            "dimension": self.dimension,
            "allowed_mask": self.allowed_mask.tolist(),
            "orthonormal_basis": self.basis.tolist(),
            "basis_note": "Basis vectors are not individual physical responses; allowed does not mean nonzero.",
            "tolerance": self.tolerance,
        }


def _solve(blocks: Sequence[FloatArray], shape: tuple[int, ...], tol: float) -> BornChargeSpace:
    columns = int(np.prod(shape))
    if not blocks:
        vectors = np.eye(columns)
    else:
        equations = np.vstack(blocks)
        # Avoid allocating a quadratic-in-row-count U for a long constraint stack.
        _, singular, vh = np.linalg.svd(equations, full_matrices=equations.shape[0] < columns)
        threshold = tol * max(1.0, float(singular[0]) if len(singular) else 1.0)
        rank = int(np.count_nonzero(singular > threshold))
        vectors = vh[rank:]
    return BornChargeSpace(shape, vectors.reshape((-1,) + shape), tol)


def atomic_born_space(operations: Iterable[Any], *, tol: float = DEFAULT_TOL) -> BornChargeSpace:
    """Solve Z=R Z R^T for ONE site, using that site's stabilizer.

    Supplying a whole-crystal group is valid only if every supplied
    operation returns the selected site to itself modulo translation.
    """
    tol = _tolerance(tol)
    matrices = _rotations(operations, tol)
    return _solve([np.eye(9) - np.kron(r, r) for r in matrices], (3, 3), tol)


def _permutations(value: Sequence[Sequence[int]], count: int) -> NDArray[np.int64]:
    array = np.asarray(value)
    if array.ndim != 2 or array.shape[0] != count or array.shape[1] < 1:
        raise ValueError("permutations must have shape (number_of_operations, number_of_atoms)")
    if not np.issubdtype(array.dtype, np.integer):
        raise ValueError("permutations must use integer, zero-based atom indices")
    if not np.all(np.sort(array, axis=1) == np.arange(array.shape[1])):
        raise ValueError("each atom mapping must be a permutation of 0..N-1")
    return array.astype(np.int64, copy=False)


def structure_born_space(
    operations: Iterable[Any],
    permutations: Sequence[Sequence[int]],
    *,
    acoustic_sum_rule: bool = True,
    tol: float = DEFAULT_TOL,
) -> BornChargeSpace:
    """Solve Z[pi_g(k)]=R_g Z[k] R_g^T for a specified atom mapping.

    The optional ASR enforces sum_k Z[k]=0 for a neutral insulating cell.
    pi_g(k) is the destination of source atom k.  Translation parts are
    encoded in pi_g and must not be discarded before computing mappings.
    This dense solver is intended for small primitive cells, not large DWs.
    """
    tol = _tolerance(tol)
    matrices = _rotations(operations, tol)
    mappings = _permutations(permutations, len(matrices))
    n = mappings.shape[1]
    blocks: list[FloatArray] = []
    for r, permutation in zip(matrices, mappings):
        block = np.zeros((9*n, 9*n))
        action = np.kron(r, r)
        for k, destination in enumerate(permutation):
            block[9*k:9*k+9, 9*destination:9*destination+9] += np.eye(9)
            block[9*k:9*k+9, 9*k:9*k+9] -= action
        blocks.append(block)
    if acoustic_sum_rule:
        blocks.append(np.tile(np.eye(9), (1, n)))
    return _solve(blocks, (n, 3, 3), tol)


def mode_born_space(
    operations: Iterable[Any],
    mode_representations: Iterable[ArrayLike],
    *,
    tol: float = DEFAULT_TOL,
) -> BornChargeSpace:
    """Solve R_g C = C D_g for delta p=e C delta q, delta q'=D_g delta q.

    D may be nonorthogonal if the mode basis is not normalized, but must
    be invertible.  The caller supplies representations of the SAME
    operations on a symmetry-closed mode subspace, in matching order.
    """
    tol = _tolerance(tol)
    matrices = _rotations(operations, tol)
    modes = [_real_array(d, "mode representation") for d in mode_representations]
    if len(modes) != len(matrices) or not modes:
        raise ValueError("one mode representation is required for each operation")
    m = modes[0].shape[0] if modes[0].ndim == 2 else 0
    if m < 1 or any(d.shape != (m, m) for d in modes):
        raise ValueError("mode representations must be square and of equal positive size")
    for d in modes:
        singular = np.linalg.svd(d, compute_uv=False)
        if singular[-1] <= tol * max(1.0, float(singular[0])):
            raise ValueError("mode representations must be invertible")
    blocks = [np.kron(r, np.eye(m)) - np.kron(np.eye(3), d.T)
              for r, d in zip(matrices, modes)]
    return _solve(blocks, (3, m), tol)


def _layer_actions(
    operations: Iterable[Any], layer_exchange: Sequence[int] | None, tol: float,
) -> tuple[list[FloatArray], list[FloatArray], list[int], list[int]]:
    matrices = _rotations(operations, tol)
    eta: list[int] = []
    for r in matrices:
        if not (np.allclose(r[2, :2], 0, atol=tol, rtol=0)
                and np.allclose(r[:2, 2], 0, atol=tol, rtol=0)
                and np.isclose(abs(r[2, 2]), 1, atol=tol, rtol=0)):
            raise ValueError("operation does not preserve the xy layer plane; "
                             "rotate into a layer-normal-z Cartesian frame first")
        eta.append(1 if r[2, 2] > 0 else -1)
    if layer_exchange is None:
        # Valid for a free, isolated two-layer object centered about its midplane.
        exchange = list(eta)
    else:
        exchange = list(layer_exchange)
        if len(exchange) != len(matrices) or any(
            isinstance(s, (bool, np.bool_)) or not isinstance(s, (int, np.integer))
            or s not in (-1, 1) for s in exchange
        ):
            raise ValueError("layer_exchange must contain one integer +1/-1 per operation")
    actions = [chi*r[:2, :2] for chi, r in zip(exchange, matrices)]
    return matrices, actions, eta, exchange


def sliding_born_space(
    operations: Iterable[Any],
    *,
    layer_exchange: Sequence[int] | None = None,
    tol: float = DEFAULT_TOL,
) -> BornChargeSpace:
    """Allowed zeta=(d p_z/d s_x, d p_z/d s_y)/e for relative sliding.

    s=u_top-u_bottom; D_g=chi_g A_g, chi=+1 preserves layer labels,
    chi=-1 exchanges them.  General constraint: D_g.T zeta=eta_g zeta.
    When layer_exchange is omitted, chi=eta=R_zz is assumed, VALID ONLY
    for an isolated free bilayer.  Substrates/3D stacks need explicit
    mappings (and may require mode_born_space instead of this 2-mode model).
    """
    tol = _tolerance(tol)
    _, actions, eta, _ = _layer_actions(operations, layer_exchange, tol)
    return _solve([d.T - sign*np.eye(2) for d, sign in zip(actions, eta)], (2,), tol)


def bilayer_mode_born_space(
    operations: Iterable[Any],
    *,
    layer_exchange: Sequence[int] | None = None,
    tol: float = DEFAULT_TOL,
) -> BornChargeSpace:
    """Full 3x2 relative-sliding charge C, including in-plane polarization."""
    tol = _tolerance(tol)
    matrices, actions, _, _ = _layer_actions(operations, layer_exchange, tol)
    return mode_born_space(matrices, actions, tol=tol)


def _charges(value: ArrayLike) -> FloatArray:
    array = _real_array(value, "Born charges")
    if array.ndim != 3 or array.shape[1:] != (3, 3) or array.shape[0] < 1:
        raise ValueError("Born charges must have shape (N, 3, 3), N>=1")
    return array


def project_mode_charges(charges: ArrayLike, displacement_patterns: ArrayLike) -> FloatArray:
    """C[i,a]=sum_(k,j) Z[k,i,j] U[k,j,a], for u=U q.

    U must include any mass normalization explicitly.  No implicit mass
    factors or rigid-layer assumptions are applied.
    """
    z = _charges(charges)
    u = _real_array(displacement_patterns, "displacement patterns")
    if u.ndim != 3 or u.shape[:2] != (len(z), 3) or u.shape[2] < 1:
        raise ValueError("displacement patterns must have shape (N, 3, number_of_modes)")
    return np.einsum("kij,kja->ia", z, u)


def bilayer_mode_charge(charges: ArrayLike, layer_labels: Sequence[int]) -> FloatArray:
    """Return C (3x2) for u_top=+s/2, u_bottom=-s/2.

    labels are +1 (top) or -1 (bottom), with both present.  This formula
    never silently imposes neutrality.  C[2] is zeta.  The simpler upper
    layer sum equals C only when the relevant ASR residual is zero.
    """
    z = _charges(charges)
    labels = np.asarray(layer_labels)
    if labels.shape != (len(z),) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("layer_labels must be an integer vector of length N")
    if set(labels.tolist()) != {-1, 1}:
        raise ValueError("both +1 (top) and -1 (bottom) labels must be present")
    u = np.zeros((len(z), 3, 2))
    u[:, 0, 0] = labels/2
    u[:, 1, 1] = labels/2
    return project_mode_charges(z, u)


def acoustic_sum_residual(charges: ArrayLike) -> FloatArray:
    """Return sum_k Z[k] without altering the data."""
    return _charges(charges).sum(axis=0)


def field_induced_forces(charges: ArrayLike, field_v_per_angstrom: ArrayLike) -> FloatArray:
    """Linear field-induced forces in eV/angstrom for Z in units of |e|.

    field input is in V/angstrom, NOT V/nm.  These are the induced force
    corrections, not the zero-field mechanical force or a switching test.
    """
    z = _charges(charges)
    field = _real_array(field_v_per_angstrom, "field")
    if field.shape != (3,):
        raise ValueError("field must have three Cartesian components")
    return np.einsum("kij,i->kj", z, field)


def read_outcar_born(path: str | Path) -> FloatArray:
    """Read the LAST raw VASP OUTCAR Born-charge block, without transposing.

    Reject an incomplete last block rather than falling back silently to
    an earlier calculation.  Check NIONS when available.  py4vasp's
    force-first convention is different and must be transposed upstream.
    """
    text = Path(path).read_text(encoding="utf-8")
    headers = list(re.finditer(r"^\s*BORN EFFECTIVE CHARGES[^\n]*$", text, re.M))
    if not headers:
        raise ValueError("OUTCAR contains no BORN EFFECTIVE CHARGES block")
    lines = text[headers[-1].end():].splitlines()
    ions: list[list[list[float]]] = []
    position = 0
    while position < len(lines):
        line = lines[position].strip()
        ion_match = re.fullmatch(r"ion\s+(\d+)", line, flags=re.I)
        if ion_match:
            if int(ion_match.group(1)) != len(ions)+1:
                raise ValueError("Born block has nonconsecutive ion indices")
            rows: list[list[float]] = []
            for expected_row in (1, 2, 3):
                position += 1
                if position >= len(lines):
                    raise ValueError("incomplete final Born-charge block")
                fields = lines[position].split()
                if len(fields) != 4 or fields[0] != str(expected_row):
                    raise ValueError("malformed or incomplete Born-charge matrix row")
                try:
                    rows.append([float(x.replace("D", "E").replace("d", "e")) for x in fields[1:]])
                except ValueError as exc:
                    raise ValueError("invalid numeric Born-charge matrix row") from exc
            ions.append(rows)
        elif line and not set(line).issubset({"-"}):
            if ions:
                break
            raise ValueError("unexpected text before first ion in Born-charge block")
        position += 1
    if not ions:
        raise ValueError("empty final Born-charge block")
    nions = re.findall(r"\bNIONS\s*=\s*(\d+)", text[:headers[-1].start()])
    if nions and len(ions) != int(nions[-1]):
        raise ValueError("final Born-charge block does not contain NIONS atoms")
    return _charges(ions)


def _named_operations(name: str | int) -> tuple[Any, ...]:
    # Import lazily: matrix-only use and tests do not load large registries.
    from .point_groups import iter_crystallographic_point_groups, point_group_operations

    if isinstance(name, str):
        key = name.strip().replace(" ", "").replace("−", "-").lower()
        # Prefer HM symbols to numeric IDs: string '3' means C3, int 3 means C2.
        matches = [g for g in iter_crystallographic_point_groups()
                   if key in {g.hm_symbol.lower(), g.schoenflies_symbol.lower(),
                              *(a.lower() for a in g.aliases)}]
        if len(matches) == 1:
            name = matches[0].number
        elif key.isdecimal():
            name = int(key)
    return point_group_operations(name)


def point_group_born_space(name: str | int, *, tol: float = DEFAULT_TOL) -> BornChargeSpace:
    """Use the existing 32-point-group registry as a SITE-group template.

    Registry settings, not point-group names alone, specify axis directions.
    For a noncanonical orientation, pass rotated Cartesian matrices to
    atomic_born_space instead.
    """
    return atomic_born_space(_named_operations(name), tol=tol)


def main(argv: Sequence[str] | None = None) -> int:
    """Independent module CLI; leaves the existing group-ops CLI unchanged."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for kind in ("atomic", "sliding"):
        child = sub.add_parser(kind)
        selection = child.add_mutually_exclusive_group(required=True)
        selection.add_argument("--point-group", help="HM or Schoenflies symbol in the toolkit setting")
        selection.add_argument("--operations", type=Path, help="JSON list of Cartesian 3x3 matrices")
        child.add_argument("--tol", type=float, default=DEFAULT_TOL)
        if kind == "sliding":
            child.add_argument("--exchange", type=Path, help="JSON list: +1 keeps, -1 exchanges layers")
    project = sub.add_parser("project")
    project.add_argument("--outcar", type=Path, required=True)
    project.add_argument("--layers", type=Path, required=True, help="JSON integer vector: +1 top, -1 bottom")
    args = parser.parse_args(argv)
    try:
        if args.command == "project":
            charges = read_outcar_born(args.outcar)
            labels = json.loads(args.layers.read_text(encoding="utf-8"))
            mode = bilayer_mode_charge(charges, labels)
            output = {"convention": CONVENTION, "mode_charge_3x2": mode.tolist(),
                      "zeta": mode[2].tolist(),
                      "acoustic_sum_residual": acoustic_sum_residual(charges).tolist()}
        else:
            operations = (_named_operations(args.point_group) if args.point_group is not None
                          else json.loads(args.operations.read_text(encoding="utf-8")))
            if args.command == "atomic":
                space = atomic_born_space(operations, tol=args.tol)
            else:
                exchange = (json.loads(args.exchange.read_text(encoding="utf-8"))
                            if args.exchange is not None else None)
                space = sliding_born_space(operations, layer_exchange=exchange, tol=args.tol)
            output = space.to_dict()
            if args.command == "atomic":
                output["scope"] = "one site, using its site stabilizer (not blindly the crystal point group)"
            else:
                output["scope"] = "direct linear Ez coupling to an xy relative-sliding mode, not switchability"
                output["layer_mapping"] = ("explicit chi" if args.exchange else
                                            "assumed chi=Rzz: isolated free bilayer only")
        print(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False))
    except (ValueError, OSError, ImportError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
