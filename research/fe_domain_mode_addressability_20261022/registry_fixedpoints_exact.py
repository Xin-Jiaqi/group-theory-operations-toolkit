"""Exact rational fixed-point solver for 2D stacking registry-space automorphisms.

Given a lattice-coordinate integer rotation R with det(I-R)!=0,
solve (I-R)x = lattice integer modulo 1, with fractions, not floats.
The count is |det(I-R)|. Fixed point existence is not enough for a
physical response zero: the operator representation must forbid it,
and the actual bilayer must realize the operation (including its Seitz
translation and layer permutation).

Example: C3 fixes 3 points, C4 fixes 2, layer-swap u->-u fixes 4.
The symmetry-forced zero-set unions are hex 6 and square 4 for the
specific layer-odd Sz/current tensor; not a novel group theorem.
No external dependencies. Usage: python registry_fixedpoints_exact.py
"""
from fractions import Fraction


def fixed_points(R, nrange=4):
    A = [[int(i == j)-R[i][j] for j in range(2)] for i in range(2)]
    d = A[0][0]*A[1][1]-A[0][1]*A[1][0]
    if not d:
        raise ValueError("Fixed locus is not isolated; use a line solver")
    result = set()
    for i in range(-nrange, nrange+1):
        for j in range(-nrange, nrange+1):
            x = Fraction(A[1][1]*i-A[0][1]*j, d) % 1
            y = Fraction(-A[1][0]*i+A[0][0]*j, d) % 1
            result.add((x, y))
    assert len(result) == abs(d)
    return sorted(result)


def gate():
    C3 = [[0,-1],[1,-1]]
    C4 = [[0,-1],[1,0]]
    layer_swap = [[-1,0],[0,-1]]
    c3 = set(fixed_points(C3))
    c4 = set(fixed_points(C4))
    minus = set(fixed_points(layer_swap))
    assert (len(c3), len(c4), len(minus)) == (3,2,4)
    assert c3 & minus == {(Fraction(0), Fraction(0))}
    assert c4 <= minus
    assert len(c3 | minus) == 6
    assert len(c4 | minus) == 4
    print("PASS C3 fixed:", sorted(c3))
    print("PASS C4 fixed:", sorted(c4))
    print("PASS layer-swap fixed:", sorted(minus))
    print("PASS hex union 6; square union 4; counts are crystallographic")


if __name__ == "__main__":
    gate()
