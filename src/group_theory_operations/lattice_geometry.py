"""Cartesian point operations of the supplied two-periodic translation lattice."""
import itertools
import numpy as np


def lattice_point_operations(lattice, *, metric_rtol=1e-6):
    """Enumerate U^T G U = G on a reduced primitive plane, including +/- normal.

    Returns Cartesian matrices, not catalog labels. No hidden standard orientation
    is assumed. Nonprimitive supercells describe their OWN translation lattice;
    callers must recover a primitive lattice to request additional crystal symmetry.
    Requires the optional materials-structure-core geometry API.
    """
    from materials_structure_core import reduced_plane_basis
    if not np.isfinite(metric_rtol) or not 0 < metric_rtol < .01:
        raise ValueError("metric_rtol must be finite, positive and below .01")
    basis, _ = reduced_plane_basis(lattice)
    A = basis.T
    G = A.T @ A
    smin = np.linalg.svd(A, compute_uv=False)[-1]
    columns = []
    for length in np.diag(G):
        bound = int(np.ceil(np.sqrt(length*(1+metric_rtol))/smin))
        if bound > 128:
            raise ValueError("lattice operation search exceeds certified budget")
        cols = []
        for pair in itertools.product(range(-bound,bound+1), repeat=2):
            v = np.array(pair)
            if np.isclose(v@G@v,length,rtol=metric_rtol,atol=1e-10): cols.append(v)
        columns.append(cols)
    n = np.cross(*basis); n /= np.linalg.norm(n)
    inverse = np.linalg.pinv(A)
    result = []
    for u,v in itertools.product(*columns):
        U = np.column_stack((u,v))
        if abs(round(np.linalg.det(U))) != 1: continue
        if not np.allclose(U.T@G@U,G,rtol=metric_rtol,atol=metric_rtol*np.min(np.diag(G))): continue
        for sign in (1,-1):
            R = A@U@inverse + sign*np.outer(n,n)
            if not np.allclose(R.T@R,np.eye(3),atol=5*metric_rtol): continue
            result.append(R)
    return tuple(sorted(result,key=lambda x: tuple(np.round(x,10).ravel())))
