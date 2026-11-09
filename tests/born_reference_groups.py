"""Independent reference generators for tests, NOT a runtime catalog.

Cartesian embedding: principal axis z, monoclinic mirror xy, secondary
proper twofold x, vertical mirror xz.  Cubic C3 cycles x->y->z->x.
The six symbolic invariant templates below are handwritten analytic
expectations, not inferred from the solver under test.
"""
from __future__ import annotations
import numpy as np

I = np.eye(3)
INV = -I
C2Z = np.diag([-1., -1., 1.])
C2X = np.diag([1., -1., -1.])
MH = np.diag([1., 1., -1.])
MV = np.diag([1., -1., 1.])
C4 = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
C3 = np.array([[-.5, -np.sqrt(3)/2, 0.], [np.sqrt(3)/2, -.5, 0.], [0., 0., 1.]])
C6 = np.array([[.5, -np.sqrt(3)/2, 0.], [np.sqrt(3)/2, .5, 0.], [0., 0., 1.]])
C3XYZ = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])

# number, HM, Schoenflies, generators, order, form
GROUPS = [
 (1, '1', 'C1', [I], 1, 'A'),
 (2, '-1', 'Ci', [INV], 2, 'A'),
 (3, '2', 'C2', [C2Z], 2, 'B'),
 (4, 'm', 'Cs', [MH], 2, 'B'),
 (5, '2/m', 'C2h', [C2Z, INV], 4, 'B'),
 (6, '222', 'D2', [C2Z, C2X], 4, 'C'),
 (7, 'mm2', 'C2v', [C2Z, MV], 4, 'C'),
 (8, 'mmm', 'D2h', [C2Z, C2X, INV], 8, 'C'),
 (9, '4', 'C4', [C4], 4, 'D'),
 (10, '-4', 'S4', [-C4], 4, 'D'),
 (11, '4/m', 'C4h', [C4, INV], 8, 'D'),
 (12, '422', 'D4', [C4, C2X], 8, 'E'),
 (13, '4mm', 'C4v', [C4, MV], 8, 'E'),
 (14, '-42m', 'D2d', [-C4, C2X], 8, 'E'),
 (15, '4/mmm', 'D4h', [C4, C2X, INV], 16, 'E'),
 (16, '3', 'C3', [C3], 3, 'D'),
 (17, '-3', 'C3i', [C3, INV], 6, 'D'),
 (18, '32', 'D3', [C3, C2X], 6, 'E'),
 (19, '3m', 'C3v', [C3, MV], 6, 'E'),
 (20, '-3m', 'D3d', [C3, C2X, INV], 12, 'E'),
 (21, '6', 'C6', [C6], 6, 'D'),
 (22, '-6', 'C3h', [-C6], 6, 'D'),
 (23, '6/m', 'C6h', [C6, INV], 12, 'D'),
 (24, '622', 'D6', [C6, C2X], 12, 'E'),
 (25, '6mm', 'C6v', [C6, MV], 12, 'E'),
 (26, '-6m2', 'D3h', [C3, C2X, MH], 12, 'E'),
 (27, '6/mmm', 'D6h', [C6, C2X, INV], 24, 'E'),
 (28, '23', 'T', [C3XYZ, C2X], 12, 'F'),
 (29, 'm-3', 'Th', [C3XYZ, C2X, INV], 24, 'F'),
 (30, '432', 'O', [C3XYZ, C4], 24, 'F'),
 (31, '-43m', 'Td', [C3XYZ, -C4], 24, 'F'),
 (32, 'm-3m', 'Oh', [C3XYZ, C4, INV], 48, 'F'),
]
FORMS = {
 'A': [['a','b','c'],['d','e','f'],['g','h','i']],
 'B': [['a','b','0'],['c','d','0'],['0','0','e']],
 'C': [['a','0','0'],['0','b','0'],['0','0','c']],
 'D': [['a','b','0'],['-b','a','0'],['0','0','c']],
 'E': [['a','0','0'],['0','a','0'],['0','0','c']],
 'F': [['a','0','0'],['0','a','0'],['0','0','a']],
}
DIMENSIONS = {'A':9,'B':5,'C':3,'D':3,'E':2,'F':1}


def closure(generators):
    """Generate finite groups independently of the production constraint solver."""
    group = [I.copy()]
    cursor = 0
    while cursor < len(group):
        for generator in generators:
            new = generator @ group[cursor]
            if not any(np.allclose(new, old, atol=1e-9, rtol=0) for old in group):
                group.append(new)
                if len(group) > 48:
                    raise AssertionError('reference generators do not close within 48 operations')
        cursor += 1
    return group


def analytic_basis(form):
    symbols = FORMS[form]
    names = sorted({s.lstrip('-') for row in symbols for s in row if s != '0'})
    matrices = []
    for name in names:
        matrix = np.array([[1. if s==name else -1. if s=='-'+name else 0.
                            for s in row] for row in symbols])
        matrices.append(matrix / np.linalg.norm(matrix))
    return np.array(matrices)
