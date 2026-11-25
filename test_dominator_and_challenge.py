"""
Pytest verification for `upstream_dominator_and_challenge.py`, the live snapshot
`dominator_challenge_snapshot.json` (`ce1560fb`), Conductor-Aspect Quasi-GRH
exact Q-Farkas duality, and Tab 7 of `docs/index.html`.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import upstream_dominator_and_challenge as udc


SNAPSHOT_PATH = Path(__file__).resolve().parent / "dominator_challenge_snapshot.json"


def _load_snapshot() -> dict:
    return json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))


def test_shrinkage_clumping_and_agents_md_debt_census():
    snap = _load_snapshot()
    sh = snap["shrinkage_and_code_health"]

    assert snap["commit_short"] == "ce1560fb"
    assert sh["nov20_original_files"] == 175187
    assert sh["current_all_lean_files"] == 58880
    assert sh["file_reduction_pct_vs_nov20"] == 66.39
    assert sh["line_reduction_pct_vs_nov20"] == 19.81
    assert sh["clumped_files_count"] == 12024
    assert sh["clump_member_markers_count"] == 63816
    assert sh["agents_md_debt_set_option_backward_files"] == 627
    assert sh["agents_md_debt_exact_by_kernel_files"] == 3


def test_challenge_5_theorem_cone_split_zeta_is_downstream_of_dirichlet():
    """
    Verify the architectural inversion in `src/Main.lean` / `src/Challenge.lean`:
    `MainConclusion` proves `DirichletLNonvanishingRightOfHalf` across 58,857 files,
    and `riemannHypothesis` is proved downstream by specializing `N = 1, chi = 1`
    via `Claim.zeta`.
    """
    snap = _load_snapshot()
    ch = snap["challenge_5_theorem_split"]
    targets = {t["theorem"]: t for t in ch["comparator_targets"]}

    assert set(targets.keys()) == {
        "LFunction_ne_zero_of_half_lt_re",
        "riemannHypothesis",
        "riemannZeta_eq_zero_re",
        "LFunction_eq_zero_re",
        "LFunction_eq_zero_re_of_isPrimitive",
    }
    assert targets["LFunction_ne_zero_of_half_lt_re"]["reachable_src_files"] == 58857
    assert targets["riemannHypothesis"]["reachable_src_files"] == 58857
    assert "Claim.zeta" in targets["riemannHypothesis"]["proof_term_in_Main"]
    assert ch["claim_strip_total_cone_files"] == 4


def test_dominator_tree_zero_side_entrances_and_pr_islands():
    snap = _load_snapshot()
    dm = snap["dominator_pr_finder"]

    assert dm["total_reachable_nodes"] == 58858
    assert dm["max_dominator_tree_depth"] == 15
    assert dm["nontrivial_dominators_ge_2_files"] == 4413
    assert dm["pr_sized_islands_25_to_600_count"] == 155

    spine = {s["module"]: s for s in dm["macro_dominator_spine"]}
    assert spine["Development.AutomorphicArithmetic.ThetaCorrespondence.ThetaPairLifts.ThirdPlaceSeries.EveryCharacterOfLiteratureOfVWin"]["dominated_files"] == 30024
    assert spine["Development.AutomorphicArithmetic.ThetaCorrespondence.ThetaPairLifts.ThirdPlaceSeries.SiegelEisensteinSecondTermHypotheses"]["dominated_files"] == 15677

    for isl in dm["top_20_pr_ready_islands"] + dm["top_12_avt_h9_encapsulated_islands"]:
        assert isl["side_entrances_to_interior"] == 0
        assert isl["deletable_internal_files"] == isl["dominated_files"] - 1
        assert isl["dominated_files"] >= 15


def test_conductor_aspect_quasi_grh_exact_two_branch_dichotomy():
    """
    Verify exact Q-Farkas + HiGHS LP agreement across all regimes of
    `theta*(tau_k, alpha_q) = max((1 + alpha_q)/2, 2 + alpha_q - tau_k/3)`.
    """
    # 1. At alpha_q = 0, reduces to Proposition 12.1: tau_k = 18/5 -> theta* = 4/5
    r0 = udc.solve_conductor_quasi_grh(Fraction(18, 5), Fraction(0, 1))
    assert r0["theta_star_exact"] == "4/5"
    assert r0["duality_gap"] < 1e-9

    # 2. Sublinear conductor alpha_q = 1/4, twist cap tau_k = 9/2 -> theta* = 2 + 1/4 - 3/2 = 3/4 < 1
    r_sub = udc.solve_conductor_quasi_grh(Fraction(9, 2), Fraction(1, 4))
    assert r_sub["theta_star_exact"] == "3/4"
    assert r_sub["nontrivial_power_saving"] is True
    assert r_sub["farkas_duals_exact"]["table_modulus_q"] == "2/3"

    # 3. Past kink tau_k = 6 > (9 + 3/4)/2 = 39/8 at alpha_q = 1/4 -> saturates at (1 + 1/4)/2 = 5/8
    r_sat = udc.solve_conductor_quasi_grh(Fraction(6, 1), Fraction(1, 4))
    assert r_sat["theta_star_exact"] == "5/8"
    assert r_sat["farkas_duals_exact"]["table_collision"] == "1/2"
    assert r_sat["farkas_duals_exact"]["codim_span_q"] == "1/2"

    # 4. Convexity barrier alpha_q = 1 -> theta* >= 1 even for large tau_k = 10
    r_conv = udc.solve_conductor_quasi_grh(Fraction(10, 1), Fraction(1, 1))
    assert r_conv["theta_star_exact"] == "1"
    assert r_conv["nontrivial_power_saving"] is False


def test_dominator_challenge_cli_snapshot_mode(capsys):
    rc = udc.main(["--snapshot", str(SNAPSHOT_PATH)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "175,187 -> 58,880 files (-66.39%)" in out
    assert "4,413 multi-file islands" in out


def test_github_pages_dashboard_contains_pr_finder_and_conductor_grh():
    html = (Path(__file__).resolve().parent / "docs" / "index.html").read_text(encoding="utf-8")
    for required in (
        "panel-prfinder",
        "175,187",
        "4,413",
        "PoincareLelongInv",
        "SiegelEisensteinSecondTermHypotheses",
        "qgrh-svg",
    ):
        assert required in html
