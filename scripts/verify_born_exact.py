#!/usr/bin/env python3
"""Independent exact-arithmetic check of all 32 Born-charge templates.

Requires SymPy only for this optional audit, not for the library/tests.
Run from the repository root: python scripts/verify_born_exact.py
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from born_reference_groups import GROUPS, FORMS, DIMENSIONS  # noqa: E402


def main() -> int:
    try:
        import sympy as s
    except ImportError as exc:
        raise SystemExit('This optional exact audit requires SymPy (not a runtime dependency).') from exc
    results = []
    for number, hm, schoenflies, generators, order, form in GROUPS:
        rotations = [s.Matrix([[s.nsimplify(float(x), [s.sqrt(3)]) for x in row]
                              for row in r]) for r in generators]
        equations = s.Matrix.vstack(*[s.eye(9)-s.kronecker_product(r, r) for r in rotations])
        kernel = equations.nullspace()
        if len(kernel) != DIMENSIONS[form]:
            raise AssertionError(f'{hm}: exact nullity is {len(kernel)}')
        symbolic = s.Matrix([[s.sympify(x) for x in row] for row in FORMS[form]])
        for r in rotations:
            if (r*symbolic*r.T-symbolic).applyfunc(s.simplify) != s.zeros(3):
                raise AssertionError(f'{hm}: analytic tensor is not exactly invariant')
        # Parameters are independent; invariance + matching nullity proves completeness.
        parameters = sorted(symbolic.free_symbols, key=str)
        coefficients = symbolic.reshape(9, 1).jacobian(parameters)
        if len(parameters) != len(kernel) or coefficients.rank() != len(kernel):
            raise AssertionError(f'{hm}: symbolic template is not an independent complete basis')
        results.append({'number':number,'hm':hm,'schoenflies':schoenflies,
                        'order':order,'form':form,'exact_dimension':len(kernel),'passed':True})
    print(json.dumps({'method':'exact Q(sqrt(3)) generator nullspaces in SymPy',
                      'groups_checked':len(results),'all_passed':True,'results':results},indent=2))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
