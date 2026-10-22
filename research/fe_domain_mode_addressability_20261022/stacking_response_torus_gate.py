"""Mathematical test: symmetry-enforced zeros of a bilayer registry response field.

F(delta) is a real 2D covector for current -> layer-odd axial Sz, in
a HEXAGONAL homobilayer family with C3, a vertical mirror, and a
horizontal layer-exchanging mirror. Delta lives on R^2/Lambda.

F is built from a symmetry-allowed Fourier scalar f via J*grad(f):
this is a representation-theory witness, NOT microscopic Kubo or a
claim that a current-induced spin susceptibility is a potential gradient.

Exact symmetries: F(C3 delta)=C3 F(delta);
F(M delta)=-M F(delta); F(-delta)=-F(delta).
Six symmetry-fixed zeros G,A,B,C,M,N, as labeled in Ji et al.,
PRL 130, 146801 (2023), Fig.2(b). For this witness,
indices = (+1,-1,-1,-1,+1,+1) and total zero (torus Euler 0).

Only the symmetry-forced ZERO LOCI are robust to any equivariant F.
The exact vortex/saddle indices and exactly-six count are NOT general:
additional zeros may emerge and zero derivatives change index.
Uses NumPy only; run: python stacking_response_torus_gate.py
"""
import numpy as np

TWOPI = 2*np.pi
A = np.array([[1., -.5], [0., np.sqrt(3)/2]])
C3f = np.array([[0, -1], [1, -1]], int)
Mf = np.array([[-1, 1], [0, 1]], int)
C3 = A@C3f@np.linalg.inv(A)
M = A@Mf@np.linalg.inv(A)
J = np.array([[0., -1.], [1., 0.]])
STAR = np.array([[1., 0.], [0., 1.], [1., -1.]])
ZEROS = {
    "G": (0., 0.),
    "A": (.5, 0.), "B": (0., .5), "C": (.5, .5),
    "M": (1/3, 2/3), "N": (2/3, 1/3),
}

def f_gradient_hessian(frac, second=0.):
    frac = np.array(frac, float)
    grad = np.zeros(2)
    hess = np.zeros((2, 2))
    for n, amplitude in ((1, 1.), (2, second)):
        for b in STAR:
            theta = TWOPI*n*(b@frac)
            grad -= amplitude*TWOPI*n*np.sin(theta)*b
            hess -= amplitude*(TWOPI*n)**2*np.cos(theta)*np.outer(b, b)
    return grad, hess

def field(frac, second=0.):
    # Correct physical Cartesian covector, not a fractional-coordinate
    # gradient. J rotates it by 90 degrees (axial mirror parity).
    grad, _ = f_gradient_hessian(frac, second)
    return J@np.linalg.inv(A).T@grad

def jacobian(frac, second=0.):
    _, hess = f_gradient_hessian(frac, second)
    return J@np.linalg.inv(A).T@hess@np.linalg.inv(A)

def gate():
    assert np.allclose(C3.T@C3, np.eye(2))
    assert np.allclose(M.T@M, np.eye(2))
    assert np.allclose(C3@C3@C3, np.eye(2))
    assert np.allclose(M@M, np.eye(2))
    rng = np.random.default_rng(20261022)
    expected = {"G": 1, "A": -1, "B": -1, "C": -1, "M": 1, "N": 1}
    for second in (0., .05):
        for x in rng.random((50, 2)):
            assert np.linalg.norm(field(C3f@x, second)-C3@field(x, second)) < 2e-12
            assert np.linalg.norm(field(Mf@x, second)+M@field(x, second)) < 2e-12
            assert np.linalg.norm(field(-x, second)+field(x, second)) < 2e-12
            for lattice in ((1, 0), (0, 1), (1, -1)):
                assert np.linalg.norm(field(x+lattice, second)-field(x, second)) < 2e-12
        indices = {}
        for name, delta in ZEROS.items():
            assert np.linalg.norm(field(delta, second)) < 1e-12
            determinant = np.linalg.det(jacobian(delta, second))
            assert abs(determinant) > 1e-5
            indices[name] = int(np.sign(determinant))
        assert indices == expected
        assert sum(indices.values()) == 0
        print("second harmonic:", second, "indices:", indices, "sum=0")
    print("PASS C3, axial mirror, layer swap, lattice periodicity and 6 zeros")
    print("ZERO claim is general conditionally; exactly-6 and indices only witness model")
    print("NO actual electronic spin response or novelty established")

if __name__ == "__main__":
    gate()
