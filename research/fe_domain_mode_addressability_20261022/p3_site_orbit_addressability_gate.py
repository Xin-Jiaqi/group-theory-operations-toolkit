"""G2: Layer-group p3 (#65), Wyckoff-site induced spin-response test.

AFLOW plane p3 (#13) uses a1=(1,0), a2=(-1/2,sqrt(3)/2)
and C3: (x,y)->(-y,x-y); p3 layer-group #65.
Generic 3d vs 1a+1b+1c. Two layers carry synthetic collinear
AFM (+z top, -z bottom), each with three local spin-z observables.
Full site representations are derived from Seitz modulo translations.
Results are symmetry-ALLOWED dimension, not electronic Kubo values.
Neither configuration is a demonstrated FE state or switching partner.

Usage: python p3_site_orbit_addressability_gate.py
Requires numpy only. Does not modify production toolkit data or APIs.
"""
import json
import numpy as np

C3_FRAC=np.array([[0,-1],[1,-1]],int)
LATTICE=np.array([[1.,-.5],[0.,np.sqrt(3)/2]])
C3_CART=LATTICE@C3_FRAC@np.linalg.inv(LATTICE)
SPECIAL=np.array([[0.,0.],[1/3,2/3],[2/3,1/3]])
GENERIC=np.array([.137,.224])


def wrap(x):
    return np.mod(x,1.)


def orbit_three(start):
    sites=[wrap(np.asarray(start,float))]
    for _ in range(2):
        sites.append(wrap(C3_FRAC@sites[-1]))
    return np.array(sites)


def seitz_permutation(sites):
    """One site may match modulo an INTEGER lattice translation."""
    P=np.zeros((len(sites),len(sites)))
    for i,x in enumerate(sites):
        mapped=wrap(C3_FRAC@x)
        hits=[j for j,y in enumerate(sites)
              if np.linalg.norm(mapped-y-np.round(mapped-y))<1e-10]
        if len(hits)!=1:raise AssertionError((i,mapped,hits))
        P[hits[0],i]=1.
    assert np.allclose(P@P@P,np.eye(len(sites)))
    return P


def hom_basis(Q,E,tol=1e-10):
    """Real orthonormal basis for Q K=K E, row-major independent columns."""
    nq,ne=Q.shape[0],E.shape[0]
    eq=[]
    for i in range(nq):
        for j in range(ne):
            B=np.zeros((nq,ne))
            B[i,j]=1.
            eq.append((Q@B-B@E).ravel())
    M=np.array(eq).T
    _,s,vh=np.linalg.svd(M,full_matrices=True)
    rank=int(np.sum(s>tol))
    return [v.reshape(nq,ne) for v in vh[rank:]]


def dim_by_character(Q,E):
    """Character inner product for cyclic C3: (1/3) sum_g chiQ*chiE."""
    return (Q.shape[0]*E.shape[0]+np.trace(Q)*np.trace(E)
            +np.trace(Q@Q)*np.trace(E@E))/3


def audit():
    assert np.allclose(C3_CART.T@C3_CART,np.eye(2))
    assert np.allclose(np.linalg.matrix_power(C3_CART,3),np.eye(2))
    general=orbit_three(GENERIC)
    Pa=seitz_permutation(general)
    Pb=seitz_permutation(SPECIAL)
    assert np.allclose(Pa,[[0,0,1],[1,0,0],[0,1,0]])
    assert np.allclose(Pb,np.eye(3))
    QA=np.kron(np.eye(2),Pa)
    QB=np.kron(np.eye(2),Pb)
    order=np.array([1,1,1,-1,-1,-1])
    assert np.allclose(QA@order,order)
    basis_a=hom_basis(QA,C3_CART)
    basis_b=hom_basis(QB,C3_CART)
    minus_a=hom_basis(Pa,C3_CART)
    minus_b=hom_basis(Pb,C3_CART)
    assert (len(basis_a),len(basis_b))==(4,0)
    assert (len(minus_a),len(minus_b))==(2,0)
    assert abs(dim_by_character(QA,C3_CART)-4)<1e-12
    assert abs(dim_by_character(QB,C3_CART))<1e-12

    for K in basis_a:
        assert np.linalg.norm(QA@K-K@C3_CART)<1e-12
        for row in ([1]*6,[1,1,1,-1,-1,-1],
                    [1,1,1,0,0,0],[0,0,0,1,1,1]):
            assert np.linalg.norm(np.array(row)@K)<1e-12
    for K in minus_a:
        assert np.linalg.norm(Pa@K-K@C3_CART)<1e-12
        assert np.linalg.norm(np.ones(3)@K)<1e-12

    local=minus_a[0]/np.max(np.abs(minus_a[0]))
    induced=local@np.array([1.,0.])
    assert np.linalg.norm(induced)>.01
    assert abs(induced.sum())<1e-12

    # Pure relabelling is a negative control: cannot change dimension.
    V=np.array([[0,1,0],[1,0,0],[0,0,1.]])
    permuted=V@Pa@V.T
    assert len(hom_basis(permuted,C3_CART))==len(minus_a)
    assert np.linalg.norm(permuted@(V@local)-(V@local)@C3_CART)<1e-12

    return {
      "layer_group":"p3 (#65); AFLOW plane p3 (#13)",
      "C3_fractional":C3_FRAC.tolist(),
      "general_3d_frac":general.tolist(),
      "special_1a_1b_1c_frac":SPECIAL.tolist(),
      "general_permutation":Pa.astype(int).tolist(),
      "special_permutation":Pb.astype(int).tolist(),
      "ordered_six_moments":order.tolist(),
      "Hom_uniform_E_to_six_local_Sz_general":len(basis_a),
      "Hom_uniform_E_to_six_local_Sz_special":len(basis_b),
      "Hom_uniform_E_to_three_layer_odd_general":len(minus_a),
      "Hom_uniform_E_to_three_layer_odd_special":len(minus_b),
      "net_spin_and_layer_sums":"forbidden in both by C3",
      "synthetic_allowed_local_Ex_pattern":np.round(induced,8).tolist(),
      "source_to_real_FE_partner":"NOT ESTABLISHED",
      "actual_MLG_and_additional_antiunitary":"NOT VERIFIED",
      "nonzero_static_or_AC_Kubo_response":"NOT PROVEN",
      "originality":"site induced representations/hidden response have prior art"
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
    print("PASS: p3 site-orbit Hom 4 vs 0; layer-odd 2 vs 0; all sums zero")
