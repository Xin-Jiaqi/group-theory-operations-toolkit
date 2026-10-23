"""p3 (#65) site-stabilizer vs full induced-orbit response, unitary C3 ONLY.

Frobenius reciprocity: Hom_G(E,Ind_H^G spin_site)
  ~= Hom_H(Res_H E, spin_site).

Both sides computed INDEPENDENTLY (6-site induced permutation vs
single-site stabilizer nullspace). No real FE, SOC corep, Kubo or magnon.
Run: python site_frobenius_gate.py (numpy only; adjacent p3 gate needed).
"""
import numpy as np
from p3_site_orbit_addressability_gate import (
    C3_CART, C3_FRAC, GENERIC, SPECIAL, orbit_three,
    seitz_permutation, hom_basis,
)

ONE=np.eye(1)

def local_hom_dim(site,local_spin):
    """Linear invariant maps E -> local spin under site stabilizer H_r."""
    site=np.asarray(site,float)
    nspin=local_spin.shape[0]
    constraints=[]
    for n in range(3):
        frac=np.linalg.matrix_power(C3_FRAC,n)
        E=np.linalg.matrix_power(C3_CART,n)
        if np.linalg.norm(frac@site-site-
                 np.round(frac@site-site))<1e-11:
            spin=np.linalg.matrix_power(local_spin,n)
            for i in range(nspin):
                for j in range(2):
                    base=np.zeros((nspin,2))
                    base[i,j]=1.
                    constraints.append((spin@base-base@E).ravel())
    assert constraints
    s=np.linalg.svd(np.array(constraints).T,compute_uv=False)
    return nspin*2-int(sum(s>1e-10))

def gate():
    for spin,name in ((ONE,"local axial Sz"),
                      (C3_CART,"local axial transverse Sxy")):
        Pgeneral=seitz_permutation(orbit_three(GENERIC))
        Pfixed=seitz_permutation(SPECIAL)
        induced=(len(hom_basis(np.kron(Pgeneral,spin),C3_CART)),
                 len(hom_basis(np.kron(Pfixed,spin),C3_CART)))
        restricted=(local_hom_dim(GENERIC,spin),
                    sum(local_hom_dim(site,spin) for site in SPECIAL))
        assert induced==restricted
        assert restricted==((2,0) if spin.shape[0]==1 else (4,6))
        print(name,"global induced",induced,"site-stabilizer",restricted)
    print("PASS Frobenius site/induced Hom at p3; not a new group theorem")
    print("NOT TESTED: actual FE partner, magnetic antiunitary, q!=0 and magnon dynamics")

if __name__=="__main__":
    gate()
