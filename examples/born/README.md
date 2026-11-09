# Born-charge examples

Run from a complete checkout with the toolkit's existing dependencies installed:

```bash
export PYTHONPATH="$(pwd)/src${PYTHONPATH:+:$PYTHONPATH}"
python -m group_theory_operations.born_effective_charges atomic --point-group C3v
python -m group_theory_operations.born_effective_charges sliding --operations examples/born/c2x.json
python -m group_theory_operations.born_effective_charges sliding --operations examples/born/c3z.json
python -m group_theory_operations.born_effective_charges sliding --operations examples/born/mz.json
python -m group_theory_operations.born_effective_charges project --outcar examples/born/OUTCAR.synthetic --layers examples/born/four_atom_layers.json
```

Expected sliding-space dimensions: C2x=1 (along x), C3z=0, mz=2. These are the isolated-free-bilayer relative-mode rules, not the atomic-site rules.

`OUTCAR.synthetic` is intentionally synthetic: it tests the raw-OUTCAR index convention and an ASR-satisfying four-atom layer projection. It is not a VASP calculation or a material prediction. Its projected zeta is (0.05, 0) in elementary-charge units. Replace layer labels with those for the actual OUTCAR atom order.

The reader takes the last Born-charge block without transposing. Use the intended local-field convention (`LEPSILON=.TRUE.`, `LRPA=.FALSE.` for the standard case); it does not convert an excluding-local-fields result into the full response. py4vasp arrays use a different index convention and must be transposed before passing to these functions.
