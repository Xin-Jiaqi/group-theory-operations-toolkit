"""One minimal research seed: relative sliding × in-plane current -> layer-odd Sz.

Idealized generators C3z, vertical mx, horizontal mz. The latter swaps
layers. Output S_-z = Sz(top)-Sz(bottom) is AXIAL and LAYER ODD.
Drive J_parallel is an in-plane POLAR current (NOT equilibrium static E).
Relative shift U=rt-rb also transforms with a minus sign when layers swap.

Zero-shift uniform J -> S_-z is forbidden under mz. The FIRST derivative
w.r.t. U has one symmetry-allowed bilinear coupling, u_x*j_y-u_y*j_x.

Pure Seitz/representation algebra. No actual material, actual nonzero
nonequilibrium Edelstein coefficient, Kubo/relaxation protocol, or novelty.
Run: python stacking_first_seed_gate.py (numpy).
"""
import numpy as np


def gate():
    c3 = np.array([[-.5,-np.sqrt(3)/2],[np.sqrt(3)/2,-.5]])
    mx = np.diag([-1.,1.])
    mz = np.eye(2)
    # generator matrices acting on (Jx,Jy), (ux,uy), scalar S_-z
    J = [c3,mx,mz]
    U = [c3,mx,-mz]
    O = [1.,-1.,-1.]
    assert O[2] == -1 and np.allclose(J[2],np.eye(2))
    assert np.allclose(U[2],-np.eye(2))

    blocks=[]
    for j,u,o in zip(J,U,O):
        columns=[]
        for row in range(2):
            for col in range(2):
                x=np.zeros((2,2))
                x[row,col]=1.
                # Invariant scalar u.T@G@j, requiring u.T U.T G J j
                # transform in the same way as the layer-odd axial Sz.
                columns.append((o*x-u.T@x@j).ravel())
        blocks.append(np.column_stack(columns))
    matrix=np.vstack(blocks)
    _,singular,Vh=np.linalg.svd(matrix,full_matrices=True)
    rank=np.linalg.matrix_rank(matrix,tol=1e-10)
    assert rank==3 and singular[-1]<1e-12
    basis=Vh[rank:].reshape(-1,2,2)
    assert len(basis)==1
    B=basis[0]
    cross=np.array([[0,1],[-1,0]],float)/np.sqrt(2.)
    assert min(np.linalg.norm(B-cross),np.linalg.norm(B+cross))<1e-12
    for j,u,o in zip(J,U,O):
        assert np.allclose(u.T@B@j,o*B)

    # Physical layer-exchange sign matters: if one wrongly sets
    # U(mz)=+I while still using O(mz)=-1, NO first derivative survives.
    wrong_blocks=blocks.copy()
    u_wrong=np.eye(2);j=mz;o=O[2]
    columns=[]
    for row in range(2):
        for col in range(2):
            x=np.zeros((2,2));x[row,col]=1.
            columns.append((o*x-u_wrong.T@x@j).ravel())
    wrong_blocks[2]=np.column_stack(columns)
    assert np.linalg.matrix_rank(np.vstack(wrong_blocks))==4
    print("PASS: zero-shift J_parallel -> layer-odd Sz forbidden")
    print("PASS: one first-order invariant (ux*jy - uy*jx); ||residual|| < 1e-12")
    print("PASS: naive polar relative shift incorrectly predicts no derivative")
    print("Open: realistic bilayer, current relaxation protocol, Kubo, SOC, measured spin")


if __name__=="__main__":
    gate()
