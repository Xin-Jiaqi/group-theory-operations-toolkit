"""G3 cheap idea gate: response-gradient selection under registry modulation.

A uniform linear electronic response K: E -> O at an exactly mirror
symmetric registry may forbid components. Its first *adiabatic* derivative
in an EXPLICITLY specified structural displacement representation U is a
tensor Lambda: (U,E)->O. Here U is chosen as ordinary in-plane polar
coordinates (no layer exchange); changing the physical stacking
representation changes selection results. No material, Berry, DFT, phonon
resonance, nonadiabatic electron-phonon vertex or new mechanism implied.

Run: python registry_gradient_gate.py (NumPy only).
"""
import numpy as np

def allowed_indices(output_sign, field_signs, distortion_signs=None):
    if distortion_signs is None:
        return [j for j,e in enumerate(field_signs) if output_sign == e]
    return [(i,j) for i,u in enumerate(distortion_signs)
            for j,e in enumerate(field_signs) if output_sign == u*e]

def mirror_gate():
    # Spatial mirror mx = diag(-1,+1,+1).
    # Axial Sz is ODD; in-plane polar E=(Ex,Ey) transforms (-,+).
    # Assume a structural registry displacement U=(ux,uy) transforms
    # as in-plane POLAR under mx; this is not valid for every bilayer
    # layer-exchange connection, hence explicit per-operation input.
    E=(-1,1); U=(-1,1); Sz=-1
    ordinary=allowed_indices(Sz,E)
    activated=allowed_indices(Sz,E,U)
    assert ordinary==[0]                         # Ex -> Sz allowed
    assert activated==[(0,1),(1,0)]              # ux*Ey, uy*Ex

    # K_y is forbidden AT ux=0, but its derivative along odd ux
    # is allowed: Ky(ux)=sin(ux), exactly Ky(-ux)=-Ky(ux).
    eps=1e-5
    val=lambda ux: np.sin(ux)
    assert abs(val(0.))<1e-14
    assert abs((val(eps)-val(-eps))/(2*eps)-1)<1e-9
    for ux in np.linspace(-.7,.7,11):
        assert abs(val(-ux)+val(ux))<1e-14

    # Counterexample: an even K_y(ux)=cos(ux) would violate covariance.
    assert abs(np.cos(.2)+np.cos(-.2))>1

    print("PASS mirror mx: static E->Sz allowed Ex only")
    print("PASS chosen POLAR sliding U: first derivatives ux*Ey and uy*Ex")
    print("PASS symmetry-fixed Ky(0)=0 and odd finite-registry activation")
    print("NOT a new effect or dynamical spin-phonon coefficient; site/layer action REQUIRED")

if __name__=="__main__":
    mirror_gate()
