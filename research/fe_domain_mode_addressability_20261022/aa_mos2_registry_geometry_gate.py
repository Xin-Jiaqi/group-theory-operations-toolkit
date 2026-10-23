#!/usr/bin/env python3
"""Geometry-only m_z covariance gate for aligned AA bilayer MoS2.

Source geometry: Venkateswarlu et al., PRB 102, 081103(R) (2020):
a=3.18 A; S height = +/-0.49115*a relative to Mo; AA Mo-Mo spacing=6.80 A.
Checks interlayer pair distances only; no SK orbital matrices, SOC, band
structure, response or material novelty claim.

For relative shift u, place the top layer at +u/2 and bottom at -u/2.
Horizontal mirror swaps layers and exchanges the upper/lower S basis inside
each 1H layer. For top basis atom a, bottom atom b and lattice vector R, the
mapped pair at -u is (top mirror(b), bottom mirror(a), -R).

Run: python aa_mos2_registry_geometry_gate.py
"""
import itertools
import numpy as np

A = 3.18
Z_S = 0.49115 * A
D_MO_MO = 6.80
A1 = np.array([A*np.sqrt(3)/2, -A/2])
A2 = np.array([A*np.sqrt(3)/2, A/2])
BASIS = [
    ("Mo", np.array([0.0, 0.0]), 0.0, 0),
    ("S", np.array([0.0, A/np.sqrt(3)]), +Z_S, 2),
    ("S", np.array([0.0, A/np.sqrt(3)]), -Z_S, 1),
]
MIRROR = {0: 0, 1: 2, 2: 1}


def pos(layer, idx, R, u):
    _, xy, z, _ = BASIS[idx]
    shift = (+0.5 if layer == "T" else -0.5) * np.asarray(u)
    zc = (+0.5 if layer == "T" else -0.5) * D_MO_MO
    return np.array([*(xy + R[0]*A1 + R[1]*A2 + shift), zc + z])


def distance(a, b, R, u):
    return np.linalg.norm(pos("T", a, (0, 0), u) - pos("B", b, R, u))


def gate():
    ulist = [np.array([0.0, 0.0]), np.array([0.037*A, 0.011*A]),
             np.array([-0.037*A, -0.011*A])]
    shifts = list(itertools.product(range(-2, 3), repeat=2))
    max_error = 0.0
    checked = 0
    for u in ulist[1:]:
        for a, b, R in itertools.product(range(3), range(3), shifts):
            d1 = distance(a, b, R, u)
            d2 = distance(MIRROR[b], MIRROR[a], (-R[0], -R[1]), -u)
            max_error = max(max_error, abs(d1-d2))
            checked += 1
    assert max_error < 1e-12, max_error
    for a in range(3):
        for b in range(3):
            assert BASIS[a][0] == BASIS[MIRROR[a]][0]
            assert BASIS[b][0] == BASIS[MIRROR[b]][0]
    d_ss = D_MO_MO - 2*Z_S
    assert abs(d_ss - 3.676) < 0.01, d_ss
    print("PASS: AA MoS2 layer-swapping m_z pair-distance covariance")
    print(f"PASS: {checked} pair/translation comparisons; max |delta d|={max_error:.3e} A")
    print(f"PASS: interfacial S-S plane separation={d_ss:.4f} A (source table ~3.676 A)")
    print("LIMIT: geometry only; no SK orbital matrices, SOC, Kubo/Boltzmann or response computed")


if __name__ == "__main__":
    gate()
