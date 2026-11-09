#!/usr/bin/env python3
"""Regenerate the canonical reference table from independent test fixtures.

These fixtures have explicit Cartesian generators; this does not infer
site symmetry from a structure.  Use check_born_registry.py to additionally
compare a full installation's registry and named-group adapter.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'tests'))
from group_theory_operations.born_effective_charges import atomic_born_space, sliding_born_space
from born_reference_groups import GROUPS, FORMS, closure, analytic_basis


def main() -> int:
    records=[]
    for number,hm,sf,generators,order,form in GROUPS:
        space=atomic_born_space(generators)
        reference=analytic_basis(form).reshape((-1,9))
        np.testing.assert_allclose(space.projector,reference.T@reference,atol=1e-9)
        full=closure(generators)
        character=sum(np.trace(r)**2 for r in full)/len(full)
        if not np.isclose(character,space.dimension):
            raise AssertionError('character check failed')
        layer_compatible=all(np.allclose(r[2,:2],0) and np.allclose(r[:2,2],0) for r in generators)
        slide=sliding_born_space(generators) if layer_compatible else None
        records.append({
            'number':number,'hm_symbol':hm,'schoenflies_symbol':sf,
            'group_order':order,'atomic_form':form,'symbolic_tensor':FORMS[form],
            'atomic_dimension':space.dimension,'atomic_allowed_mask':space.allowed_mask.tolist(),
            'cartesian_generators':[r.tolist() for r in generators],
            'layer_plane_preserved_in_this_setting':layer_compatible,
            'isolated_bilayer_sliding_dimension':None if slide is None else slide.dimension,
            'isolated_bilayer_sliding_projector':None if slide is None else slide.projector.tolist(),
            'dimension_from_group_character':round(character),
        })
    data={
        'schema_version':1,
        'convention':'row=polarization/electric field; column=displacement/force',
        'setting':'principal axis z; monoclinic mirror xy; secondary twofold x; explicit generators authoritative',
        'atomic_scope':'ONE atomic site and its stabilizer, not the total neutral-cell charge tensor',
        'sliding_scope':'free isolated bilayer in xy, relative mode s=u_top-u_bottom, chi=eta=Rzz',
        'warnings':['A point-group name without orientation is not sufficient for a real material.',
                    'Nonzero allowed entries are not predicted values.',
                    'Null sliding response at one configuration is not a global switching prohibition.'],
        'forms':FORMS,'point_groups':records,
    }
    output=ROOT/'docs'/'born-charge-reference.json'
    output.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(output)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
