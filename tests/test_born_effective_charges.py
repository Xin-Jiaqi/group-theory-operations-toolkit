"""Regression tests for static Born-charge selection rules and conventions."""
from __future__ import annotations
import contextlib
import importlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose

from group_theory_operations import born_effective_charges as bec
from born_reference_groups import (
    GROUPS, DIMENSIONS, I, INV, C2Z, C2X, MH, MV, C3, C4, C6,
    closure, analytic_basis,
)


class TestAll32PointGroups(unittest.TestCase):
    pass


def _make_group_tests(record):
    number, hm, schoenflies, generators, expected_order, form = record
    def nullspace(self):
        answer = bec.atomic_born_space(generators)
        self.assertEqual(answer.dimension, DIMENSIONS[form])
        reference = analytic_basis(form).reshape((-1, 9))
        assert_allclose(answer.projector, reference.T @ reference, atol=1e-9)
        for tensor in answer.basis:
            for r in generators:
                assert_allclose(tensor, r @ tensor @ r.T, atol=1e-9)
    def group_average(self):
        group = closure(generators)
        self.assertEqual(len(group), expected_order)
        reynolds = sum(np.kron(r, r) for r in group)/len(group)
        assert_allclose(reynolds @ reynolds, reynolds, atol=1e-9)
        assert_allclose(reynolds.T, reynolds, atol=1e-9)
        answer = bec.atomic_born_space(generators)
        assert_allclose(reynolds, answer.projector, atol=1e-9)
        assert_allclose(bec.atomic_born_space(group).projector, answer.projector, atol=1e-9)
        dimension_from_character = sum(np.trace(r)**2 for r in group)/len(group)
        self.assertAlmostEqual(dimension_from_character, DIMENSIONS[form], places=9)
    def rotation_covariance(self):
        rng = np.random.default_rng(2100+number)
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        rotated = bec.atomic_born_space([q @ r @ q.T for r in generators])
        original = bec.atomic_born_space(generators)
        action = np.kron(q, q)
        assert_allclose(rotated.projector, action @ original.projector @ action.T, atol=1e-9)
    def bilayer_classification(self):
        # Only these canonical embeddings preserve the xy plane; cubic is rejected.
        if number > 27:
            with self.assertRaises(ValueError):
                bec.sliding_born_space(generators)
            return
        result = bec.sliding_born_space(generators)
        expected = 2 if hm in {'1', 'm'} else 0
        self.assertEqual(result.dimension, expected)
        full = bec.bilayer_mode_born_space(generators)
        # The zeta block of the full invariant projector is the same subspace.
        assert_allclose(full.projector[4:6, 4:6], result.projector, atol=1e-9)
    return nullspace, group_average, rotation_covariance, bilayer_classification


for record in GROUPS:
    for name, function in zip(
        ('atomic_nullspace', 'group_average_and_character', 'rotation_covariance', 'bilayer_classification'),
        _make_group_tests(record),
    ):
        setattr(TestAll32PointGroups, f'test_{record[0]:02d}_{name}', function)


class TestScientificEdgeCases(unittest.TestCase):
    def test_c3_allows_antisymmetric_xy(self):
        tensor = np.array([[0., 2., 0.],[-2., 0., 0.],[0., 0., 0.]])
        assert_allclose(bec.atomic_born_space([C3]).project(tensor), tensor, atol=1e-9)
        assert_allclose(bec.atomic_born_space([C3, MV]).project(tensor), 0, atol=1e-9)

    def test_mirror_site_vs_layer_exchange(self):
        site = bec.atomic_born_space([MH])
        self.assertFalse(site.allowed_mask[2, :2].any())
        self.assertEqual(bec.sliding_born_space([MH]).dimension, 2)
        structure = bec.structure_born_space([MH], [[1, 0]])
        self.assertEqual(structure.dimension, 4)  # mixed z/in-plane entries survive ASR.
        for z in structure.basis:
            assert_allclose(z[1], MH @ z[0] @ MH, atol=1e-9)
            assert_allclose(z.sum(axis=0), 0, atol=1e-9)
        projected = np.array([bec.bilayer_mode_charge(z, [1, -1])[2] for z in structure.basis])
        self.assertEqual(np.linalg.matrix_rank(projected, tol=1e-9), 2)

    def test_inversion_site_vs_layer_exchange(self):
        self.assertEqual(bec.atomic_born_space([INV]).dimension, 9)
        self.assertEqual(bec.sliding_born_space([INV]).dimension, 0)
        self.assertEqual(bec.structure_born_space([INV], [[1, 0]]).dimension, 0)
        self.assertEqual(bec.structure_born_space([INV], [[1, 0]], acoustic_sum_rule=False).dimension, 9)

    def test_neutral_single_atom_asr(self):
        self.assertEqual(bec.structure_born_space([I], [[0]]).dimension, 0)
        self.assertEqual(bec.atomic_born_space([I]).dimension, 9)

    def test_c3_general_site_orbit_not_site_stabilizer(self):
        raw = np.arange(9, dtype=float).reshape(3, 3)
        orbit = np.array([np.linalg.matrix_power(C3, i) @ raw @ np.linalg.matrix_power(C3.T, i)
                          for i in range(3)])
        answer = bec.structure_born_space([C3], [[1, 2, 0]], acoustic_sum_rule=False)
        self.assertEqual(answer.dimension, 9)
        assert_allclose(answer.project(orbit), orbit, atol=1e-8)
        self.assertNotEqual(orbit[0,2,0], 0)
        assert_allclose(orbit.sum(axis=0)[2,:2], 0, atol=1e-8)
        assert_allclose(bec.atomic_born_space([C3]).project(raw)[2,:2], 0, atol=1e-8)

    def test_atomic_nonzero_can_cancel_within_layer(self):
        z = np.zeros((4,3,3))
        z[:,2,0] = [1., -1., 2., -2.]
        assert_allclose(bec.bilayer_mode_charge(z, [1,1,-1,-1]), 0)

    def test_top_only_formula_requires_asr(self):
        z = np.arange(36, dtype=float).reshape(4,3,3)
        c = bec.bilayer_mode_charge(z, [1,1,-1,-1])
        upper = z[:2,:,:2].sum(axis=0)
        assert_allclose(upper-c, .5*bec.acoustic_sum_residual(z)[:,:2])
        z[-1] -= z.sum(axis=0)
        assert_allclose(bec.bilayer_mode_charge(z,[1,1,-1,-1]), z[:2,:,:2].sum(axis=0))

    def test_force_index_and_units(self):
        z = np.zeros((1,3,3)); z[0,2,0] = .02; z[0,0,2] = 7.
        # 0.3 V/nm = 0.03 V/angstrom.
        assert_allclose(bec.field_induced_forces(z,[0,0,.03]), [[.0006,0,0]])
        assert_allclose(bec.field_induced_forces(z,[.03,0,0]), [[0,0,.21]])

    def test_mode_projection_virtual_work(self):
        rng = np.random.default_rng(10)
        z = rng.normal(size=(5,3,3)); u=rng.normal(size=(5,3,4)); e=rng.normal(size=3)
        c=bec.project_mode_charges(z,u)
        forces=bec.field_induced_forces(z,e)
        assert_allclose(c.T @ e, np.einsum('kj,kja->a',forces,u))

    def test_diagonal_atomic_bec_can_drive_nonrigid_mode(self):
        z=np.array([np.diag([1.,2.,3.]),np.diag([-1.,-2.,-3.])])
        u=np.zeros((2,3,1)); u[:,2,0]=[.5,-.5]
        self.assertEqual(bec.project_mode_charges(z,u)[2,0],3.)

    def test_in_plane_c2_selects_line(self):
        assert_allclose(bec.sliding_born_space([C2X]).projector, np.diag([1.,0.]),atol=1e-9)
        assert_allclose(bec.sliding_born_space([MV]).projector, np.diag([1.,0.]),atol=1e-9)
        assert_allclose(bec.sliding_born_space([C2X,MH]).projector, np.diag([1.,0.]),atol=1e-9)

    def test_two_distinct_mirror_axes_remove_line(self):
        q=np.array([[2**-.5,-2**-.5,0],[2**-.5,2**-.5,0],[0,0,1.]])
        self.assertEqual(bec.sliding_born_space([MV,q@MV@q.T]).dimension,0)

    def test_rotated_allowed_sliding_direction(self):
        angle=.381; q=np.array([[np.cos(angle),-np.sin(angle),0],
                               [np.sin(angle),np.cos(angle),0],[0,0,1]])
        direction=q[:2,0]
        result=bec.sliding_born_space([q@C2X@q.T])
        assert_allclose(result.projector,np.outer(direction,direction),atol=1e-9)

    def test_explicit_layer_parity_changes_answer(self):
        self.assertEqual(bec.sliding_born_space([MH],layer_exchange=[-1]).dimension,2)
        self.assertEqual(bec.sliding_born_space([MH],layer_exchange=[1]).dimension,0)
        self.assertEqual(bec.sliding_born_space([INV],layer_exchange=[1]).dimension,2)

    def test_general_mode_equivariance(self):
        rng=np.random.default_rng(8)
        for r in [C3,C2X,MH,INV]:
            for chi in (-1,1):
                d=chi*r[:2,:2]
                space=bec.mode_born_space([r],[d])
                for c in space.basis:
                    assert_allclose(r@c,c@d,atol=1e-9)
                    q=rng.normal(size=2)
                    assert_allclose(r@(c@q),c@(d@q),atol=1e-9)

    def test_nonorthogonal_mode_basis(self):
        v=np.array([[2.,.4],[0.,.7]])
        d=np.linalg.inv(v)@C2X[:2,:2]@v
        result=bec.mode_born_space([C2X],[d])
        self.assertEqual(result.dimension,3)
        for c in result.basis:
            assert_allclose(C2X@c,c@d,atol=1e-9)

    def test_atom_covariance_and_asr_projection(self):
        rng=np.random.default_rng(13)
        result=bec.structure_born_space([C2X],[[1,0,3,2]])
        z=result.project(rng.normal(size=(4,3,3)))
        for k,dest in enumerate([1,0,3,2]):
            assert_allclose(z[dest],C2X@z[k]@C2X.T,atol=1e-9)
        assert_allclose(z.sum(axis=0),0,atol=1e-9)

    def test_structure_projection_matches_mode_selection(self):
        rng=np.random.default_rng(31)
        for r in [I,C2Z,C3,MH,MV,C2X,INV]:
            chi=1 if r[2,2]>0 else -1
            pi=[0,1] if chi==1 else [1,0]
            structural=bec.structure_born_space([r],[pi])
            z=structural.project(rng.normal(size=(2,3,3)))
            # For C3, the chosen fixed sites can have more symmetry than general sites;
            # this is a valid mapping model, not an arbitrary actual crystal.
            c=bec.bilayer_mode_charge(z,[1,-1])
            full=bec.bilayer_mode_born_space([r])
            assert_allclose(full.project(c),c,atol=1e-9)

    def test_endpoints_zero_does_not_imply_no_path_coupling(self):
        q=np.linspace(0.,1.,101)
        p=np.cos(np.pi*q); zeta=-np.pi*np.sin(np.pi*q)
        self.assertAlmostEqual(zeta[0],0)
        self.assertAlmostEqual(zeta[-1],0)
        self.assertGreater(np.max(np.abs(zeta)),3)
        self.assertAlmostEqual(p[-1]-p[0],-2)
        # No force at q=0; curvature of -E*p is E*pi^2 and changes with E.
        self.assertLess(1.+(-1.)*np.pi**2,0)


class TestLiteratureExamplesAndCharacters(unittest.TestCase):
    """Independent examples from Fachin et al., npj Comput. Mater. 11, 307 (2025)."""

    def test_germanene_site_3m_and_inversion_plus_neutrality(self):
        # C3v fixes each site, inversion exchanges two sites.  Each atomic
        # tensor is diag(a,a,c); equal inversion partners + ASR imply zero.
        self.assertEqual(bec.atomic_born_space([C3,MV]).dimension,2)
        with_inversion=bec.structure_born_space([C3,MV,INV],[[0,1],[0,1],[1,0]])
        self.assertEqual(with_inversion.dimension,0)
        # Breaking inversion while retaining the two local C3v sites
        # allows opposite, nonzero two-parameter tensors.
        without_inversion=bec.structure_born_space([C3,MV],[[0,1],[0,1]])
        self.assertEqual(without_inversion.dimension,2)
        for z in without_inversion.basis:
            assert_allclose(z[0],-z[1],atol=1e-9)
            assert_allclose(z[0],np.diag([z[0,0,0],z[0,0,0],z[0,2,2]]),atol=1e-9)

    def test_jacutingaite_m_and_2_over_m_axis_x(self):
        # The paper's representative sites have unique axis/normal x:
        # Zxy, Zxz, Zyx, Zzx vanish; yz and zy need not be equal.
        mx=-C2X
        for generators in ([mx],[C2X,INV]):
            result=bec.atomic_born_space(generators)
            self.assertEqual(result.dimension,5)
            self.assertEqual(result.allowed_mask.tolist(),
                             [[True,False,False],[False,True,True],[False,True,True]])
            z=np.array([[1.,0.,0.],[0.,2.,3.],[0.,4.,5.]])
            assert_allclose(result.project(z),z,atol=1e-9)

    def test_full_atom_character_and_asr_nullity(self):
        # Formula requires a FULL group with consistent permutations,
        # not merely the generator list used by the nullspace solver.
        cases=[([I,C3,C3@C3],[[0,1,2],[1,2,0],[2,0,1]]),
               ([I,C2X],[[0,1,2,3],[1,0,3,2]]),
               ([I,MH],[[0,1],[1,0]]),
               ([I,INV],[[0,1],[1,0]])]
        for group,maps in cases:
            n=len(maps[0])
            char=sum(sum(int(i==j) for i,j in enumerate(pi))*np.trace(r)**2
                     for r,pi in zip(group,maps))/len(group)
            tensor_char=sum(np.trace(r)**2 for r in group)/len(group)
            raw=bec.structure_born_space(group,maps,acoustic_sum_rule=False)
            neutral=bec.structure_born_space(group,maps)
            self.assertAlmostEqual(raw.dimension,char,places=9)
            self.assertAlmostEqual(neutral.dimension,char-tensor_char,places=9)
            for z in neutral.basis:
                assert_allclose(z.sum(axis=0),np.zeros((3,3)),atol=1e-9)


class TestValidationAndIO(unittest.TestCase):
    def test_invalid_rotations(self):
        bads=[[], [np.ones((2,2))], [np.diag([1,1,2])], [np.eye(3)*np.nan], [np.eye(3)*1j]]
        for bad in bads:
            with self.subTest(bad=str(bad)), self.assertRaises(ValueError):
                bec.atomic_born_space(bad)

    def test_fractional_hexagonal_rotation_rejected(self):
        fractional=np.array([[0,-1,0],[1,-1,0],[0,0,1]])
        with self.assertRaisesRegex(ValueError,'Cartesian'):
            bec.atomic_born_space([fractional])
        l=np.array([[1,-.5,0],[0,np.sqrt(3)/2,0],[0,0,1]])
        self.assertEqual(bec.atomic_born_space([l@fractional@np.linalg.inv(l)]).dimension,3)

    def test_bad_tolerance(self):
        for tol in [0,-1,np.nan,.1,np.inf]:
            with self.subTest(tol=tol), self.assertRaises(ValueError):
                bec.atomic_born_space([I],tol=tol)

    def test_bad_atom_mappings(self):
        for maps in [[[0,0]],[[1.,0.]],[[0,2]],[],[[True,False]]]:
            with self.subTest(maps=maps), self.assertRaises(ValueError):
                bec.structure_born_space([I],maps)

    def test_bad_modes(self):
        for modes in [[],[np.zeros((2,2))],[np.ones((2,3))]]:
            with self.subTest(modes=str(modes)), self.assertRaises(ValueError):
                bec.mode_born_space([I],modes)

    def test_bad_layer_labels(self):
        z=np.zeros((2,3,3))
        for labels in [[1,1],[1],[1.,-1.],[0,1],[True,False]]:
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                bec.bilayer_mode_charge(z,labels)

    def test_bad_layer_exchange(self):
        for exchange in [[0],[True],[1.],[],[1,-1]]:
            with self.subTest(exchange=exchange), self.assertRaises(ValueError):
                bec.sliding_born_space([I],layer_exchange=exchange)

    def test_layer_normal_mixing_rejected(self):
        tilt=np.array([[1.,0,0],[0,0,-1],[0,1,0]])
        with self.assertRaisesRegex(ValueError,'layer plane'):
            bec.sliding_born_space([tilt])

    def test_operation_record_duck_typing(self):
        op=types.SimpleNamespace(matrix_cartesian=C3)
        self.assertEqual(bec.atomic_born_space([op]).dimension,3)

    def test_zero_dim_space_json_and_projection(self):
        result=bec.sliding_born_space([C3])
        self.assertEqual(result.basis.shape,(0,2))
        self.assertEqual(result.to_dict()['allowed_mask'],[False,False])
        assert_allclose(result.project([1.,2.]),0)
        json.dumps(result.to_dict(),allow_nan=False)

    def test_basis_is_read_only(self):
        result=bec.atomic_born_space([I])
        with self.assertRaises(ValueError):
            result.basis[0,0,0]=2

    def test_outcar_parser_orientation_and_last_block(self):
        first='BORN EFFECTIVE CHARGES (including local field effects)\n---\nion 1\n1 1 2 3\n2 4 5 6\n3 7 8 9\n'
        last='BORN EFFECTIVE CHARGES (including local field effects)\n---\nion 1\n1 10 20 30\n2 40 50 60\n3 70 80 9.0D1\n'
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'OUTCAR'; path.write_text('NIONS = 1\n'+first+'other output\n'+last)
            z=bec.read_outcar_born(path)
            assert_allclose(z[0],[[10,20,30],[40,50,60],[70,80,90]])
            path.write_text('NIONS = 1\n'+first+'other output\n'+last.rsplit('\n3',1)[0])
            with self.assertRaises(ValueError):
                bec.read_outcar_born(path)

    def test_outcar_nions_and_missing_block(self):
        block='BORN EFFECTIVE CHARGES\nion 1\n1 1 2 3\n2 4 5 6\n3 7 8 9\n'
        for contents in ['NIONS=2\n'+block,'no block',block.replace('7 8 9','7 nan 9')]:
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as directory:
                path=Path(directory)/'OUTCAR';path.write_text(contents)
                with self.assertRaises(ValueError):
                    bec.read_outcar_born(path)

    def test_matrix_only_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ops.json';path.write_text(json.dumps([C3.tolist()]))
            for mode,expected in [('atomic',3),('sliding',0)]:
                stdout=io.StringIO()
                with contextlib.redirect_stdout(stdout):
                    self.assertEqual(bec.main([mode,'--operations',str(path)]),0)
                self.assertEqual(json.loads(stdout.getvalue())['dimension'],expected)

    def test_named_registry_adapter_contract(self):
        # This is an interface/contract test, not a fetched upstream-registry test.
        records=[types.SimpleNamespace(number=3,hm_symbol='2',schoenflies_symbol='C2',aliases=()),
                 types.SimpleNamespace(number=16,hm_symbol='3',schoenflies_symbol='C3',aliases=())]
        selected=[]
        stub=types.ModuleType('group_theory_operations.point_groups')
        stub.iter_crystallographic_point_groups=lambda: iter(records)
        def resolve(name):
            selected.append(name)
            return (types.SimpleNamespace(matrix_cartesian=C3 if name==16 else C2Z),)
        stub.point_group_operations=resolve
        with patch.dict(sys.modules,{'group_theory_operations.point_groups':stub}):
            self.assertEqual(bec.point_group_born_space('3').dimension,3)
            self.assertEqual(selected[-1],16)
            self.assertEqual(bec.point_group_born_space(3).dimension,5)
            self.assertEqual(selected[-1],3)
            self.assertEqual(bec.point_group_born_space('C3').dimension,3)


if __name__=='__main__':
    unittest.main()
