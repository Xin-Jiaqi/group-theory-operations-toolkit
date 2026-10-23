"""Correct transformation of relative in-plane stacking displacement under layer exchange.

An explicitly geometric, finite-dimensional Seitz-point-operation gate, no
claim of SOC/Kubo physics. Run python physical_layer_exchange_gate.py.

For a layer swapping spatial operation R, u = r_top - r_bottom transforms
as (-1)*R_parallel*u (up to registry lattice shifts), not as E_parallel.
Layer-odd spin O_- = S_top - S_bottom transforms by an extra (-1)
relative to ordinary axial spin. K0 and dK/du are classified separately.

Scope: orthogonal diagonal operations mx/mz/C2x with separately verified
layer exchange; not generic non-symmorphic ops with q phases, not a genuine
material magnetic operation, not a transport or equilibrium response.
"""
import numpy as np


def classify(R, layer_sign, spin_axis):
    """layer_sign=-1 when the operation swaps top and bottom."""
    R = np.asarray(R, dtype=int)
    assert np.array_equal(R.T @ R, np.eye(3, dtype=int))
    assert layer_sign in (-1, 1)
    E = R[:2, :2]                  # in-plane polar external E
    U = layer_sign * E             # relative slide r_T-r_B
    A = np.linalg.det(R) * R       # spin: axial 3-vector
    output = layer_sign * round(A[spin_axis, spin_axis])
    e = np.diag(E).astype(int).tolist()
    u = np.diag(U).astype(int).tolist()
    K0 = [j for j,ej in enumerate(e) if output == ej]
    gradient = [(i,j) for i,ui in enumerate(u)
                      for j,ej in enumerate(e) if output == ui*ej]
    return dict(E=e, u=u, layer_odd_spin=int(output),
                K0_allowed=K0, gradient_allowed=gradient)


def gate():
    mz = np.diag([1,1,-1])
    c2x = np.diag([1,-1,-1])
    a = classify(mz,-1,2)         # m_z: O=S_T^z-S_B^z
    assert a["K0_allowed"] == []
    assert a["gradient_allowed"] == [(0,0),(0,1),(1,0),(1,1)]
    b = classify(mz,-1,0)         # m_z: O=S_T^x-S_B^x
    assert b["K0_allowed"] == [0,1]
    assert b["gradient_allowed"] == []
    c = classify(c2x,-1,2)
    assert c["K0_allowed"] == [0]
    assert c["gradient_allowed"] == [(0,1),(1,0)]

    # WRONG shortcut: replace U by the external polar-field matrix for m_z.
    # Keep the actual layer-odd output sign fixed. This would falsely
    # classify every derivative of O_-^z as forbidden.
    wrong_u = np.diag(mz[:2,:2]).astype(int).tolist()
    wrong_gradient = [(i,j) for i,ui in enumerate(wrong_u)
                            for j,ej in enumerate(a["E"])
                            if a["layer_odd_spin"] == ui*ej]
    assert wrong_u == [1,1]
    assert wrong_gradient == []

    K = lambda ux,uy: np.array([[.2*ux+.3*uy, -.8*ux+1.1*uy]])
    for ux,uy in [(0.,0.),(.2,-.03),(-.8,.17)]:
        assert np.allclose(K(-ux,-uy),-K(ux,uy))
    print("PASS m_z layer-odd Sz",a)
    print("PASS m_z layer-odd Sx",b)
    print("PASS C2x layer-odd Sz",c)
    print("NAIVE U AS POLAR IS WRONG",wrong_gradient)
    print("No inference of T-breaking, FE stability or nonequilibrium Kubo response")


if __name__ == "__main__":
    gate()
