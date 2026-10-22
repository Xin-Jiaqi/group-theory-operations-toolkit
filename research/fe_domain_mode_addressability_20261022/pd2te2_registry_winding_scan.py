#!/usr/bin/env python3
"""Parameter/radius scan of layer-even/odd current-spin response winding.

Reads the user's local BiDB SQLite through the existing read-only model.
The source database is not copied or uploaded. This is a toy-model stress
test, NOT a material prediction or gauge-covariant Kubo calculation.

Example:
  OPENBLAS_NUM_THREADS=1 python pd2te2_registry_winding_scan.py \
    --db /path/to/bidb.db --id 980 --nk 8 \
    --mu 1.0 1.5 2.0 --soc 0.24 0.48 0.72 \
    --radius 0.125 --angles 12

The integer is a sampled contour winding; always inspect min |F|, repeat
with denser k/angle meshes, and do not interpret it as robust if |F| gets
small anywhere on the contour. A finite-radius index can change when
response zeros cross the loop.
"""
import argparse
import json
import numpy as np
from pd2te2_te_p_soc_boltzmann_gate import Model


def calc_winding(model, radius, mu, soc, nk, nangles, T):
    fields = [[], []]
    for i in range(nangles):
        theta = 2*np.pi*i/nangles
        u = radius*np.array([np.cos(theta), np.sin(theta)])
        chi = model.chi_current(u, nk=nk, mu=mu, T=T, soc=soc)
        for channel in (0, 1):
            fields[channel].append(chi[channel, 0] + 1j*chi[channel, 1])
    out = {}
    for channel, name in ((0, "layer_odd"), (1, "layer_even")):
        z = np.asarray(fields[channel], dtype=complex)
        # A winding is undefined if the sampled field vanishes; do not
        # divide by zero or report a spurious integer in that case.
        norms = np.abs(z)
        sampled_zero = bool(np.any(norms < 1e-10))
        if sampled_zero:
            winding = None
        else:
            increments = np.angle(np.roll(z, -1)/z)
            winding = float(np.sum(increments)/(2*np.pi))
        out[name] = {
            "winding": winding,
            "min_norm": float(np.min(norms)),
            "max_norm": float(np.max(norms)),
            "sampled_zero": sampled_zero,
        }
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--id", type=int, choices=[980, 975], default=980)
    p.add_argument("--nk", type=int, default=8)
    p.add_argument("--mu", type=float, nargs="+", default=[1.0, 1.5, 2.0])
    p.add_argument("--soc", type=float, nargs="+", default=[0.24, 0.48, 0.72])
    p.add_argument("--radius", type=float, nargs="+", default=[0.125])
    p.add_argument("--angles", type=int, default=12)
    p.add_argument("--T", type=float, default=0.13)
    a = p.parse_args()
    model = Model(a.db, a.id)
    for radius in a.radius:
        for mu in a.mu:
            for soc in a.soc:
                result = calc_winding(model, radius, mu, soc,
                                      a.nk, a.angles, a.T)
                print(json.dumps({
                    "bidb_id": a.id, "nk": a.nk,
                    "angles": a.angles, "radius_angstrom": radius,
                    "mu_eV": mu, "soc_eV": soc,
                    "T_eV": a.T, **result
                }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
