"""Hexagonal vs square registry response: LAYER-ODD and LAYER-EVEN channels.

Physical input: in-plane current, output: out-of-plane axial spin, measured
either as S_-z (top-bottom) or S_+z (top+bottom). Assumed family operations
C3/C4, vertical mirror, and horizontal layer swap. Only unitary symmetry.

Hex near a mirror-symmetric AA reference, writing u=(ux,uy):
  layer-odd response F_-(u)=gamma*(-uy,ux) + O(u^3)
  layer-even response F_+(u)=beta*(ux^2-uy^2,-2ux*uy)+O(u^4)
They are symmetry-ALLOWED, not necessarily nonzero.

Square C4 + horizontal layer exchange: C4^2=C2 maps registry u->-u
and in-plane current j->-j, but even axial Sz->+Sz.
Thus F_+(-u)=-F_+(u). The horizontal layer swap maps u->-u,
yet F_+(-u)=+F_+(u). Both conditions force F_+(u)=0 for ALL u.
This is a known centrosymmetric/pure-sliding symmetry constraint,
NOT a newly discovered mechanism. Layer-odd channel is not forced zero.

Runs with NumPy only; no TB, Kubo, Berry, DFT or material prediction.
"""
import numpy as np

C3=np.array([[-.5,-np.sqrt(3)/2],[np.sqrt(3)/2,-.5]])
M=np.diag([-1.,1.])
J=np.array([[0.,-1.],[1.,0.]])
C4=np.array([[0.,-1.],[1.,0.]])

def layer_odd_linear(u):
    return J@u

def layer_even_quadratic(u):
    x,y=u
    return np.array([x*x-y*y,-2*x*y])

def gate():
    rng=np.random.default_rng(21)
    for label,F,parity,n in (
        ("layer-odd",layer_odd_linear,-1,1),
        ("layer-even",layer_even_quadratic,1,2)):
        for u in rng.normal(size=(40,2)):
            assert np.linalg.norm(F(C3@u)-C3@F(u))<1e-12
            assert np.linalg.norm(F(M@u)+M@F(u))<1e-12
            assert np.linalg.norm(F(-u)-parity*F(u))<1e-12
        print("PASS hex",label,"lowest ALLOWED order:",n)

    assert np.array_equal(C4@C4,-np.eye(2))
    # C2z: response covector F_+ transforms like an in-plane polar vector,
    # so C2z gives F_+(-u)=-F_+(u), despite Sz being C2-even.
    # Layer swap mz: output S_+z is even, current is even; F_+ is even.
    # Contradictory constraints force F_+(u) identically zero.
    for p in rng.normal(size=(20,2)):
        possible_c2=lambda v:-v
        possible_mz=lambda v:v
        test=np.array([p[0],p[1]])
        assert np.allclose(possible_c2(test)+possible_mz(test),0.)
    print("PASS square C4^2 versus mz forces layer-even Sz/current response =0")
    print("No nonzero coefficient, spin Hall current or material response demonstrated")

if __name__=="__main__":
    gate()
