"""Antiunitary-domain Onsager test for finite-frequency spectral Kubo.

Independent two/four-band experiments, no real material, no new theorem.
Run: python antiunitary_kubo_onsager_gate.py (numpy only)

e^{-i omega t}, retarded omega+En-Em+i*eta, eta>0; full eigenspaces.
"""
import numpy as np

I=np.eye(2,dtype=complex)
X=np.array([[0.,1.],[1.,0.]],complex)
Y=np.array([[0.,-1j],[1j,0.]],complex)
Z=np.diag([1.,-1.]).astype(complex)


def chi(H,A,B,omega=.6,eta=.03,nocc=1):
    """Occupied/nonoccupied spectral Lehmann sum, no guessed DC protocol."""
    e,u=np.linalg.eigh(H)
    a=u.conj().T@A@u
    b=u.conj().T@B@u
    f=np.arange(len(e))<nocc
    de=e[:,None]-e[None,:]
    df=f[:,None].astype(float)-f[None,:].astype(float)
    return np.sum(df*a*b.T/(omega+de+1j*eta))


def antiunitary_map(u,h):
    return u@h.conj()@u.conj().T


def spinless():
    h=.5*X+.8*Y+.7*Z
    a,b=Z,X
    ca=chi(h,a,b)
    swapped=chi(h,b,a)
    partner=chi(h.conj(),a,b)
    assert abs(partner-swapped)<1e-13
    assert abs(partner-ca)>.1
    assert abs(partner-np.conj(ca))>.01
    print("K-only: source chi_AB",ca,"source chi_BA",swapped,
          "mapped-domain chi_AB",partner)


def spinful():
    # True spin-1/2 time reversal, theta^2=-1; no TR symmetry assumed for H.
    u=np.kron(I,1j*Y)
    assert np.linalg.norm(u@u.conj()+np.eye(4))<1e-13
    rng=np.random.default_rng(712)
    r=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    h=(r+r.conj().T)/2
    hb=antiunitary_map(u,h)
    a=np.kron(I,Z)/2       # physical spin z, eta_A=-1
    b=np.kron(X,I)         # orbital polarization, eta_B=+1
    assert np.linalg.norm(antiunitary_map(u,a)+a)<1e-13
    assert np.linalg.norm(antiunitary_map(u,b)-b)<1e-13
    for w in (0.,.24,.8):
        direct=chi(h,a,b,omega=w,nocc=2)
        recip=chi(h,b,a,omega=w,nocc=2)
        mapped=chi(hb,a,b,omega=w,nocc=2)
        assert abs(mapped+recip)<1e-12
        if w>0: assert abs(mapped+direct)>1e-4
        print("Theta^2=-1, omega",w,"source spin<-P",direct,
              "reciprocal P<-spin",recip,"mapped spin<-P",mapped)


if __name__=="__main__":
    spinless()
    spinful()
    print("PASS: chi_AB^partner(omega) = eta_A eta_B chi_BA^source(omega)")
    print("FAIL for naive same-index sign and naive complex conjugation")
    print("Scope: finite-model Lehmann only; magnetic tensor coreps still OPEN")
