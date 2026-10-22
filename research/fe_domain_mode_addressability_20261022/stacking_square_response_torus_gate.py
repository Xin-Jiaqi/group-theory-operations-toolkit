"""Second family: square homobilayer registry response (C4 + mirror + layer swap).

A real current-to-layer-odd-axial-Sz response covector F(delta) is
constrained by F(C4*d)=C4*F(d), F(mx*d)=-mx*F(d), F(-d)=-F(d).
Symmetry-forced zeros at all four half-translation registries; C4 fixed
registries form a subset. Lowest-harmonic Fourier witness F=J*grad(Phi),
Phi=cos(2pi x)+cos(2pi y), has indices (+,-,-,+).
Not a material/Kubo/new-theorem claim.

Run: python stacking_square_response_torus_gate.py (numpy).
"""
import numpy as np

T=2*np.pi
J=np.array([[0.,-1.],[1.,0.]])
C4=np.array([[0.,-1.],[1.,0.]])
M=np.diag([-1.,1.])
FIXED={"G":(0.,0.),"A":(.5,0.),"B":(0.,.5),"C":(.5,.5)}

def field(frac, harmonic=0.):
    p=np.asarray(frac,float)
    grad=-T*np.sin(T*p)-harmonic*2*T*np.sin(2*T*p)
    return J@grad

def jac(frac, harmonic=0.):
    p=np.asarray(frac,float)
    h=-T*T*np.cos(T*p)-harmonic*(2*T)**2*np.cos(2*T*p)
    return J@np.diag(h)

def gate():
    rng=np.random.default_rng(20261022)
    expected={"G":1,"A":-1,"B":-1,"C":1}
    for harmonic in (0.,.05):
        for p in rng.random((30,2)):
            assert np.linalg.norm(field(C4@p,harmonic)-C4@field(p,harmonic))<1e-12
            assert np.linalg.norm(field(M@p,harmonic)+M@field(p,harmonic))<1e-12
            assert np.linalg.norm(field(-p,harmonic)+field(p,harmonic))<1e-12
            assert np.linalg.norm(field(p+np.array([1,0]),harmonic)-field(p,harmonic))<1e-12
        indices={}
        for name,p in FIXED.items():
            assert np.linalg.norm(field(p,harmonic))<1e-12
            det=np.linalg.det(jac(p,harmonic))
            assert abs(det)>1e-5
            indices[name]=int(np.sign(det))
        assert indices==expected
        assert sum(indices.values())==0
        print("square harmonic",harmonic,"indices",indices,"sum 0")
    print("PASS square four forced zeros, C4 and mirror covariance, index budget")
    print("No spin Hamiltonian, material response, exact-zero-count theorem or claim of novelty")

if __name__=="__main__":
    gate()
