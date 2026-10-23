"""Symmetry-only layer-EVEN axial Sz/current response in hexagonal registry space.

In the same C3, vertical mirror and layer-swap family as the layer-odd
response, the layer-EVEN observable S_+z is EVEN under horizontal
layer swap: F_+(-delta)=+F_+(delta), while vertical mirror still gives
F_+(M delta)=-M F_+(delta) and C3 vector covariance.

At AA/G: C3 forbids a constant; layer-swap evenness forbids linear
terms. Generic first allowed term is quadratic, complex-coordinate
F_x+i F_y proportional to (u_x-i u_y)^2, with winding -2.
The three C3 fixed points G,M,N must be zero; the A/B/C half lattice
registries are not forced zeros by layer exchange for this observable.

Witness F=J grad(sin(2pi x)-sin(2pi y)-sin(2pi(x-y)), plus optional
higher harmonic. Pure representation theory, NOT a material, Kubo or
a newly discovered vortex theorem.

Run: python registry_layer_even_spin_gate.py (numpy).
"""
import numpy as np

T=2*np.pi
A=np.array([[1.,-.5],[0.,np.sqrt(3)/2]])
B=np.linalg.inv(A)
C3F=np.array([[0,-1],[1,-1]],int)
MF=np.array([[-1,1],[0,1]],int)
C3=A@C3F@B
M=A@MF@B
J=np.array([[0.,-1.],[1.,0.]])
STAR=np.array([[1.,0.],[0.,-1.],[-1.,1.]])
FORCED={"G":(0.,0.),"M":(1/3,2/3),"N":(2/3,1/3)}

def field_and_jac(frac, second=0.):
    p=np.asarray(frac, float)
    grad=np.zeros(2)
    hess=np.zeros((2,2))
    for n,coeff in ((1,1.),(2,second)):
        for k in STAR:
            theta=T*n*(k@p)
            grad += coeff*T*n*np.cos(theta)*k
            hess -= coeff*(T*n)**2*np.sin(theta)*np.outer(k,k)
    return J@B.T@grad, J@B.T@hess@B

def gate():
    rng=np.random.default_rng(44)
    for second in (0., .05):
        for p in rng.random((25,2)):
            f=field_and_jac(p,second)[0]
            assert np.linalg.norm(field_and_jac(C3F@p,second)[0]-C3@f)<1e-12
            assert np.linalg.norm(field_and_jac(MF@p,second)[0]+M@f)<1e-12
            assert np.linalg.norm(field_and_jac(-p,second)[0]-f)<1e-12
        for p in FORCED.values():
            assert np.linalg.norm(field_and_jac(p,second)[0])<1e-12
        for name in ("M","N"):
            assert np.linalg.det(field_and_jac(FORCED[name],second)[1])>100
        theta=np.linspace(0,2*np.pi,1001)
        radius=.01
        circle=np.array([
            field_and_jac(B@np.array([radius*np.cos(t),radius*np.sin(t)]),second)[0]
            for t in theta])
        assert np.min(np.linalg.norm(circle,axis=1))>1e-4
        angles=np.unwrap(np.arctan2(circle[:,1],circle[:,0]))
        winding=(angles[-1]-angles[0])/(2*np.pi)
        assert abs(winding+2)<1e-9
        print("PASS layer-even second harmonic",second,
              "index(G)=-2, index(M)=index(N)=+1, total=0")
    print("This is a synthetic response field, not computed physical spin accumulation.")

if __name__=="__main__":
    gate()
