import unittest
import numpy as np
from group_theory_operations.lattice_geometry import lattice_point_operations
from group_theory_operations.stacking import polarization_space


class LatticeGeometryTests(unittest.TestCase):
    def test_counts_and_rebasis(self):
        for base,count in [(np.array([[3.,0,0],[0,4,0]]),8),
                           (np.array([[3.,0,0],[0,3,0]]),16),
                           (np.array([[3.,0,0],[1.5,3*np.sqrt(3)/2,0]]),24),
                           (np.array([[3.,0,0],[.7,4.,0]]),4),
                           (np.array([[2.,1.,0],[-2.,1.,0]]),8)]:
            expected=lattice_point_operations(base)
            self.assertEqual(len(expected),count)
            U=np.array([[1,9],[0,1]])
            changed=lattice_point_operations(U@base)
            self.assertEqual(len(changed),count)
            for op in expected:
                self.assertTrue(any(np.allclose(op,x,atol=1e-8) for x in changed))

    def test_rotation_and_polarity(self):
        q,_=np.linalg.qr(np.random.default_rng(5).normal(size=(3,3)))
        base=np.array([[3.,0,0],[0,4,0]])
        original=lattice_point_operations(base)
        rotated=lattice_point_operations(base@q.T)
        for op in original:
            self.assertTrue(any(np.allclose(q@op@q.T,x) for x in rotated))
        p=polarization_space([q@np.diag([-1.,-1.,1.])@q.T],normal=q[:,2])
        self.assertEqual(p.polar_type,'OP')
