#!/usr/bin/env python3
"""Check the additive module against the installed toolkit's actual registry.

This complements the independent local tests.  It requires a complete
checkout/install, not merely the additive patch directory.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'tests'))
import numpy as np
from group_theory_operations.born_effective_charges import point_group_born_space
from group_theory_operations.point_groups import iter_crystallographic_point_groups
from born_reference_groups import GROUPS, analytic_basis


def main() -> int:
    records=list(iter_crystallographic_point_groups())
    expected={record[0]:record for record in GROUPS}
    if len(records)!=32:
        raise AssertionError('the toolkit registry must contain 32 point groups')
    rows=[]
    for group in records:
        number,hm,_,_,_,form=expected[group.number]
        result=point_group_born_space(number)
        by_name=point_group_born_space(group.hm_symbol)
        basis=analytic_basis(form).reshape((-1,9))
        np.testing.assert_allclose(result.projector,basis.T@basis,atol=1e-8)
        np.testing.assert_allclose(by_name.projector,result.projector,atol=1e-8)
        rows.append({'number':number,'hm':hm,'dimension':result.dimension,'passed':True})
    print(json.dumps({'upstream_registry_groups':len(rows),'all_passed':True,'rows':rows},indent=2))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
