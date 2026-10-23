"""Falsification-first B+D idea gate. Python standard library ONLY.

This is a conservative CLAIM gate, not a scientific novelty oracle or a
group-theory theorem prover. User-authored input cannot independently certify
whether a source or symmetry mapping is valid. This script NEVER runs VASP,
touches confidential source files, submits jobs, or changes GitHub status.

Usage: python idea_elimination_gate.py IDEA.json
       python idea_elimination_gate.py --selftest
"""
import argparse
import json
import sys

KINDS = {"unitary_full", "antiunitary_full", "structural_only",
         "nonconjugate_magnetic", "not_applicable", "unknown"}
MECHANISMS = {"new_intrinsic_channel", "same_forward_ac_sign",
              "hidden_local_response", "moire_shift_current_vortex",
              "mode_projected_switching", "geometry_texture", "other"}
REQUIRED = ("defined_operators", "response_protocol_defined",
            "source_symmetry_checked", "physical_switch_pair_verified",
            "magnetic_site_map_checked", "novelty_equation_compared",
            "experimentally_distinguishable")
TESTS = ("unitary_covariance", "antiunitary_onsager_swap",
         "site_relabel_invariance", "rank_vs_hom", "real_vs_complex_ac",
         "toy_off_controls", "gauge_and_contact", "peer_same_observable")

# Specific 2025-26 overlap reports. Not comprehensive; never assume anything
# without the operator-level comparison specified in the actual record.
KNOWN_OVERLAPS = {
    "10.1103/dmzg-ck2t": "80-layer-group FE magnetic connection",
    "10.1103/qhvs-pfkq": "layer Edelstein",
    "10.1103/np2g-86h3": "FE charge-spin switching",
    "10.1038/s41699-025-00552-6": "moire shift-photocurrent vortices",
    "10.1038/s41467-025-60647-y": "experimental moire polarization merons",
    "10.1103/PhysRevB.98.245129": "magnetic multipoles",
    "10.1103/PhysRevB.102.165143": "hidden local physical effects",
    "10.1103/ktpw-7s2b": "AFM SOT symmetry classification",
    "10.1038/s41467-026-69635-2": "FE-spin symmetry classes",
}


def evaluate(d):
    """Return reproducible verdict with ALL failures, never automatically pass novelty."""
    problems = []
    eliminations = []
    if not isinstance(d, dict):
        raise ValueError("Input must be a JSON object")
    for k in ("idea_id", "hypothesis", "observable", "domain", "claim",
              "evidence", "closest_peers", "controls"):
        if k not in d:
            raise ValueError("Missing field: " + k)
    obs, domain, claim, evidence = (
        d["observable"], d["domain"], d["claim"], d["evidence"])
    if not all(isinstance(x, dict) for x in (obs, domain, claim, evidence)):
        raise ValueError("observable, domain, claim, evidence must be objects")
    kind = domain.get("connection_type", "unknown")
    mech = claim.get("mechanism")
    if kind not in KINDS or mech not in MECHANISMS:
        raise ValueError("Unknown connection_type or mechanism")
    if not all(obs.get(k) for k in ("source", "measurement", "spatial_support",
                                     "response_order", "time_protocol")):
        problems.append("BLOCKED_OPERATOR_DEFINITION")
    if kind == "unknown":
        problems.append("BLOCKED_REAL_CONNECTION")
    if kind == "unitary_full" and mech == "new_intrinsic_channel":
        eliminations.append(
            "EXACT_UNITARY_COVARIANCE: mapped full-Hamiltonian domains "
            "cannot gain independent channels by coordinate rotation")
    if kind == "antiunitary_full" and mech == "same_forward_ac_sign":
        eliminations.append(
            "ANTIUNITARY_ONSAGER: retarded reciprocity exchanges SOURCE "
            "and MEASUREMENT; same-forward sign claim is unsupported")
    if mech == "hidden_local_response" and claim.get("generic_novelty", False):
        eliminations.append(
            "HIDDEN_SITE_PRIOR_ART: induced Wyckoff response and compensated "
            "local magnetic multipoles have established predecessors")
    if mech == "moire_shift_current_vortex" and claim.get("generic_novelty", False):
        eliminations.append(
            "SHIFT_CURRENT_VORTEX_PRIOR_ART: npj 2D Mater. Appl. 9, 43 "
            "(2025), doi:10.1038/s41699-025-00552-6")
    if claim.get("mode_count_equals_hom_dimension", False):
        eliminations.append(
            "HOM_NOT_MODE_COUNT: allowed coefficient space dimension is "
            "not physical normal-mode count or rank of a given response")
    if claim.get("assumes_full_FE_partner_changes_wyckoff_multiplicity", False):
        eliminations.append(
            "WYCKOFF_CONJUGACY: exact symmetry partner preserves orbit "
            "stabilizer conjugacy and orbit multiplicity")
    peers = d["closest_peers"]
    if not isinstance(peers, list):
        raise ValueError("closest_peers must be a list")
    for p in peers:
        if not isinstance(p, dict) or not p.get("doi"):
            raise ValueError("Each peer needs a DOI")
        if p.get("same_observable") and p.get("same_mechanism") and \
           p.get("relation") == "subsumed":
            eliminations.append("PRIOR_ART_SUBSUMED:" + p["doi"])
        if p.get("relation") not in ("subsumed", "different", "unresolved"):
            raise ValueError("peer relation must be subsumed/different/unresolved")
    for key in REQUIRED:
        if evidence.get(key) is not True:
            problems.append("MISSING_EVIDENCE:" + key)
    controls = d["controls"]
    if not isinstance(controls, dict):
        raise ValueError("controls must be object")
    for key in TESTS:
        if controls.get(key) not in ("pass", "fail", "not_applicable",
                                     "not_tested"):
            raise ValueError("control " + key + ": expected pass/fail/not_applicable/not_tested")
    if kind == "unitary_full" and controls["unitary_covariance"] != "pass":
        problems.append("UNITARY_DOMAIN_COVARIANCE_NOT_VERIFIED")
    if kind == "antiunitary_full" and \
       controls["antiunitary_onsager_swap"] != "pass":
        problems.append("RECIPROCAL_KUBO_NOT_VERIFIED")
    if controls["site_relabel_invariance"] == "fail":
        problems.append("GAUGE_OR_SITE_RELABEL_FAILURE")
    if controls["rank_vs_hom"] == "fail":
        problems.append("MISTAKEN_MODE_COUNT")
    if claim.get("material_prediction"):
        if evidence.get("native_numeric_convergence") is not True:
            problems.append("NO_ACCEPTED_MATERIAL_NUMERICS")
        if evidence.get("real_response_vertices") is not True:
            problems.append("NO_SAME_GAUGE_RESPONSE_VERTICES")
    if not evidence.get("closest_peer_equations_verified", False):
        problems.append("CLOSEST_PEER_EQUATION_REVIEW_PENDING")
    if any(p.get("relation") == "unresolved" for p in peers):
        problems.append("UNRESOLVED_CLOSE_PEER")
    if not peers:
        problems.append("NO_CLOSEST_PEER")
    status = "NO_GO" if eliminations else (
        "BLOCKED" if problems else "PHYSICS_CANDIDATE")
    return {"idea_id": d["idea_id"], "status": status,
            "hard_eliminations": sorted(set(eliminations)),
            "missing_or_failed_gates": sorted(set(problems)),
            "not_a_novelty_or_material_certificate": True,
            "new_compute_authorized": False,
            "known_overlap_coverage": {
                p["doi"]: KNOWN_OVERLAPS[p["doi"]]
                for p in peers if p["doi"] in KNOWN_OVERLAPS}}


def fixture(mechanism, kind="structural_only"):
    return {"idea_id": "synthetic", "hypothesis": "A test only",
            "observable": {"source": "E_x", "measurement": "site S_z",
                           "spatial_support": "sites",
                           "response_order": "linear",
                           "time_protocol": "static"},
            "domain": {"connection_type": kind},
            "claim": {"mechanism": mechanism, "generic_novelty": False},
            "evidence": {},
            "closest_peers": [{"doi": "10.1103/qhvs-pfkq",
                               "relation": "unresolved",
                               "same_observable": False,
                               "same_mechanism": False}],
            "controls": {k: "not_tested" for k in TESTS}}


def selftest():
    x = fixture("new_intrinsic_channel", "unitary_full")
    assert evaluate(x)["status"] == "NO_GO"
    x = fixture("same_forward_ac_sign", "antiunitary_full")
    assert evaluate(x)["status"] == "NO_GO"
    x = fixture("hidden_local_response")
    x["claim"]["generic_novelty"] = True
    assert evaluate(x)["status"] == "NO_GO"
    x = fixture("moire_shift_current_vortex")
    x["claim"]["generic_novelty"] = True
    assert evaluate(x)["status"] == "NO_GO"
    x = fixture("mode_projected_switching")
    assert evaluate(x)["status"] == "BLOCKED"
    x["claim"]["mode_count_equals_hom_dimension"] = True
    assert evaluate(x)["status"] == "NO_GO"
    x["claim"]["mode_count_equals_hom_dimension"] = False
    x["claim"]["material_prediction"] = True
    r = evaluate(x)
    assert r["status"] == "BLOCKED"
    assert "NO_ACCEPTED_MATERIAL_NUMERICS" in r["missing_or_failed_gates"]
    x = fixture("mode_projected_switching", "not_applicable")
    x["evidence"] = {k: True for k in REQUIRED}
    x["evidence"]["closest_peer_equations_verified"] = True
    x["closest_peers"][0]["relation"] = "different"
    x["controls"] = {k: "not_applicable" for k in TESTS}
    assert evaluate(x)["status"] == "PHYSICS_CANDIDATE"
    assert evaluate(x)["not_a_novelty_or_material_certificate"]
    print("PASS 8 falsification gates, BLOCKED, NO_GO and noncertifying CANDIDATE")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="?")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
    elif args.input:
        with open(args.input, encoding="utf-8") as f:
            print(json.dumps(evaluate(json.load(f)), indent=2, ensure_ascii=False))
    else:
        parser.error("Supply idea JSON or --selftest")
