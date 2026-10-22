"""Test UNITARY symmetry-domain response covariance and glide phases.

Toy algebra only: not complete 80/528-layer-group enumeration, not an
antiunitary finite-frequency Kubo statement, not a material prediction.
Run: python domain_covariance_gate.py  (numpy)
"""
import numpy as np


def equivariant_dimension(mode_ops, field_ops, tol=1e-10):
    """Real dimension of D_mode(g) K = K D_field(g)."""
    nq, ne = mode_ops[0].shape[0], field_ops[0].shape[0]
    equations = []
    for q, e in zip(mode_ops, field_ops):
        for i in range(nq):
            for j in range(ne):
                unit = np.zeros((nq, ne))
                unit[i, j] = 1.
                equations.append((q@unit-unit@e).ravel())
    return nq*ne - int(np.linalg.matrix_rank(np.array(equations).T, tol=tol))


def transform(K, U_mode, U_field):
    """Only if the two FULL magnetic states are EXACT unitary partners."""
    return U_mode@K@U_field.T


def gate():
    # Mirror mx forbids Ex but allows Ey to excite an even scalar mode.
    E_mirror=np.diag([-1.,1.])
    Q_even=np.eye(1)
    assert equivariant_dimension([Q_even],[E_mirror]) == 1
    assert equivariant_dimension([Q_even],[np.eye(2)]) == 2
    K=np.array([[0.,2.]])
    assert np.allclose(Q_even@K,K@E_mirror)

    # C4-related domain: lab-fixed Ex becomes bright; this is NOT an
    # intrinsic new channel, since the field itself must be rotated.
    C4=np.array([[0.,-1.],[1.,0.]])
    new=transform(K,Q_even,C4)
    mirror_new=C4@E_mirror@C4.T
    assert equivariant_dimension([Q_even],[mirror_new]) == 1
    assert np.allclose(Q_even@new,new@mirror_new)
    assert np.allclose(new,[[-2.,0.]])
    assert abs(K[0,0])<1e-12 and abs(new[0,0])>1.
    assert np.allclose(new@(C4@np.array([0.,1.])),K@[0.,1.])
    assert np.allclose(np.linalg.svd(K,compute_uv=False),
                       np.linalg.svd(new,compute_uv=False))

    # Full-rank toy illustrates invariance of ALL singular values.
    rng=np.random.default_rng(2048)
    A=rng.normal(size=(2,2))
    Q=np.diag([1.,-1.])
    B=transform(A,Q,C4)
    assert np.linalg.matrix_rank(A)==np.linalg.matrix_rank(B)
    assert np.allclose(np.linalg.svd(A,compute_uv=False),
                       np.linalg.svd(B,compute_uv=False))

    # Glide along x by a/2, exchanging two sites, squares to T_a.
    # At k_x=pi/a the glide eigenvalues are +/-i, not +/-1.
    for k in (0.,np.pi/3,np.pi):
        G=np.array([[0.,np.exp(-1j*k)],[1.,0.]])
        assert np.allclose(G@G,np.exp(-1j*k)*np.eye(2))
        assert np.allclose(G.conj().T@G,np.eye(2))
        if np.isclose(k,np.pi):
            assert all(np.min(abs(np.linalg.eigvals(G)-z))<1e-12
                       for z in (1j,-1j))

    print("PASS conjugated stabilizer Hom dimension, response rank/singular values")
    print("PASS C4 lab-fixed OFF->ON is only domain orientation:", K, new)
    print("PASS glide g(k)^2=exp(-ik*a), zone-edge eigenvalues +/-i")
    print("NOT COVERED: antiunitary finite-omega, physical eigenmodes, 80/528 catalog")


if __name__=="__main__":
    gate()
