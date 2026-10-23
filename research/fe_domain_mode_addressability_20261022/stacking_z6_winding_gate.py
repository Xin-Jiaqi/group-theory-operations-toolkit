"""Minimal mathematical gate: layer-parity selection of response winding modulo six.

For a REAL (DC) planar current-to-axial-Sz response covector F(z) on the
registry disk, assume the full unperturbed bilayer family has:
  C3: F(exp(2pi i/3) z)=exp(2pi i/3) F(z);
  layer-swap mz: F(-z)= -F(z) for layer ODD, +F(z) for layer EVEN Sz;
  vertical mx: F(-conj(z))=conj(F(z)).
If F never vanishes on a C3+swap-invariant circular loop, its winding
must satisfy w_odd = 1 (mod 6), w_even = 4 (mod 6), irrespective of
coefficients. These are standard equivariant homotopy deductions, not
an original topology theorem, a material response or Kubo.

The two-leading-monomial even response beta conj(z)^2 + c z^4, with
real beta and c, can host six +1 satellites at |z|=sqrt(|beta/c|);
the origin contributes -2, so outer contour winding +4.
This mathematical witness is NOT evidence of satellite vortices in real
Pd2Te2 or its unfitted 24-band toy. The latter has mesh/cutoff sensitivity.
Run: python stacking_z6_winding_gate.py
"""
import numpy as np

ROT=np.exp(2j*np.pi/3)

def winding(f,r=.2,n=480):
    theta=2*np.pi*np.arange(n)/n
    values=f(r*np.exp(1j*theta))
    if np.any(np.abs(values)<1e-12):
        raise ValueError("Response zero lies on the integration contour")
    increments=np.angle(np.roll(values,-1)/values)
    assert np.max(np.abs(increments))<np.pi/2
    return round(float(increments.sum()/(2*np.pi)))

def symmetry_check(f,sign):
    z=np.array([.11+.19j,-.31+.26j,.2-.17j])
    assert np.allclose(f(ROT*z),ROT*f(z))
    assert np.allclose(f(-z),sign*f(z))
    assert np.allclose(f(-np.conjugate(z)),np.conjugate(f(z)))

def gate():
    odd_1=lambda z:1j*z
    odd_minus5=lambda z:1j*np.conjugate(z)**5
    even_minus2=lambda z:np.conjugate(z)**2
    even_4=lambda z:z**4
    for f,sign,w in ((odd_1,-1,1),(odd_minus5,-1,-5),
                     (even_minus2,1,-2),(even_4,1,4)):
        symmetry_check(f,sign)
        assert winding(f)==w
        assert w%6==(1 if sign==-1 else 4)

    beta=-.02
    c=1.
    f=lambda z:beta*np.conjugate(z)**2+c*z**4
    symmetry_check(f,1)
    radius=np.sqrt(abs(beta/c))
    assert winding(f,r=.09)==-2
    assert winding(f,r=.22)==4
    roots=radius*np.exp(1j*np.arange(6)*np.pi/3)
    assert np.max(np.abs(f(roots)))<1e-13
    for z in roots:
        dx=1e-5
        dfdx=(f(z+dx)-f(z-dx))/(2*dx)
        dfdy=(f(z+1j*dx)-f(z-1j*dx))/(2*dx)
        det=np.imag(np.conjugate(dfdx)*dfdy)
        assert det>0
    print("PASS: layer-odd winding = 1 mod 6; layer-even = 4 mod 6")
    print("PASS: six +1 satellite zeros, origin -2, outer contour +4")
    print("No new general theorem or microscopic material claim established")

if __name__=="__main__":
    gate()
