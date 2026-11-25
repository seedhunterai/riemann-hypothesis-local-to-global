"""
End-to-end pytest verification for the 3-Stage Open-Problem System:
  - Stage 1 (`lean_ast_ingestor.py` + `rh_blueprint.lean`): `.lean` parsing,
    `lean4export` JSON roundtrip, Skolem telescope poset & `sorry` leak catching
  - Stage 2 (`tropical_mip_auditor.py`): Disjunctive MILP hypothesis selection,
    exact rational (`Fraction`) Farkas certificates, & Guth-Maynard 30/13 saddle
  - Stage 3 (`quasi_rh_extractor.py`): Three-regime Window-Free Quasi-RH dichotomy
    and automatic slack triage of `sorry`-bearing hypotheses
"""

from __future__ import annotations

from pathlib import Path

import lean_ast_ingestor as ingest
import quasi_rh_extractor as qrh
import tropical_mip_auditor as mip


BLUEPRINT_FILE = Path(__file__).resolve().parent / "rh_blueprint.lean"


def test_stage1_lean_blueprint_and_export_json_roundtrip():
    """Verify `.lean` parsing, `lean4export` JSON roundtrip, and `sorry` uniformity checking."""
    src = BLUEPRINT_FILE.read_text(encoding="utf-8")
    decls = ingest.parse_lean_file_text(src)
    assert len(decls) == 16

    payload = ingest.export_decls_to_lean4_json(decls)
    assert payload["format"] == "lean4export-ast-v1"
    roundtrip = ingest.ingest_lean4_json(payload)
    rep = ingest.analyze_ingested_decls(roundtrip)

    assert rep["quantifier_poset"]["is_well_founded_dag"] is True
    assert rep["quantifier_poset"]["uniformity_H12_valid"] is True
    assert len(rep["extracted_asymptotic_terms"]) == 9
    assert "H9_galois_cohomology_debt" in rep["sorry_bearing_declarations"]
    assert "Lemma_1_2_window_landau" in rep["proved_declarations"]

    # Inject an illegal uniformity leak into a `sorry` theorem (`∃ (C_unif : Real)` after `∀ (x : Real)`)
    bad_src = src.replace(
        "∀ (eps : Real), ∀ (x : Real), ∃ (m_steps : Nat), True",
        "∀ (eps : Real), ∀ (x : Real), ∃ (C_unif : Real), ∃ (m_steps : Nat), True",
    )
    bad_rep = ingest.analyze_ingested_decls(ingest.parse_lean_file_text(bad_src))
    assert bad_rep["quantifier_poset"]["uniformity_H12_valid"] is False


def test_stage2_guth_maynard_zero_density_exact_saddle():
    """
    Verify exact rational primal/dual certification on both:
      - 1972 Ingham + Huxley barrier: sigma* = 3/4, A* = 12/5, gap = 7/12
      - 2024 Ingham + Huxley + Guth-Maynard: sigma* = 7/10, A* = 30/13, gap = 17/30,
        with Huxley receiving shadow price 0 at the 30/13 bottleneck!
    """
    z72 = mip.solve_zero_density_polyhedron(include_guth_maynard=False)
    assert z72["exact_rational_saddle_certified"] is True
    assert z72["worst_case_sigma_exact"] == "3/4"
    assert z72["zero_density_exponent_A_star_exact"] == "12/5"
    assert z72["prime_gap_exponent_exact"] == "7/12"

    z24 = mip.solve_zero_density_polyhedron(include_guth_maynard=True)
    assert z24["exact_rational_saddle_certified"] is True
    assert z24["worst_case_sigma_exact"] == "7/10"
    assert z24["zero_density_exponent_A_star_exact"] == "30/13"
    assert z24["prime_gap_exponent_exact"] == "17/30"
    assert z24["dual_shadow_prices_exact"]["Huxley_1972 (3/(3*sigma-1))"] == "0"
    assert z24["dual_shadow_prices_exact"]["Ingham_1940 (3/(2-sigma))"] == "1/2"
    assert z24["dual_shadow_prices_exact"]["Guth_Maynard_2024 (15/(5*sigma+3))"] == "1/2"


def test_stage2_disjunctive_milp_drops_windows_for_quasi_rh():
    """
    Verify that the MILP hypothesis selector automatically drops windowed H9 (cost 10)
    in favor of window-free H9 (cost 1) when targeting Quasi-RH (theta <= 0.85),
    while selecting windowed H9 when targeting Full RH (theta <= 0.50).
    """
    rep = qrh.run_full_open_problem_report()["stage2_disjunctive_milp_selection"]

    mq = rep["target_Quasi_RH_theta_0.85"]
    assert mq["feasible"] is True
    assert mq["requires_windows"] is False
    assert mq["selected_hypotheses"] == ["H9_window_free_tau3.6"]
    assert abs(mq["achieved_theta"] - 0.80) < 1e-5

    mf = rep["target_Full_RH_theta_0.50"]
    assert mf["feasible"] is True
    assert mf["requires_windows"] is True
    assert mf["selected_hypotheses"] == ["H9_windowed_x400"]
    assert abs(mf["achieved_theta"] - 0.50) < 1e-5


def test_stage3_three_regime_quasi_rh_dichotomy_and_slack_triage():
    """
    Verify all three regimes of the Sha-free, window-free Quasi-RH dichotomy
    with zero duality gap and exact rational Farkas certificates:
      - tau_k = 4.5 -> Regime I   (Full RH, theta* = 1/2)
      - tau_k = 3.6 -> Regime II  (Window-Free Quasi-RH, theta* = 4/5 = 0.80 < 1,
                                   and H1/H2/H4/H5/H6 all drop to SLACK_DEGRADABLE!)
      - tau_k = 2.7 -> Regime III (Farkas Impossibility theta* = 11/10 >= 1, no power saving)
    """
    r_full = qrh.extract_quasi_rh_at_tau(4.5)
    assert r_full["regime"] == "REGIME_I_FULL_WINDOW_FREE_RH"
    assert r_full["farkas_certified_in_Q"] is True
    assert r_full["farkas_exact_rational_lower_bound"] == "1/2"

    r_quasi = qrh.extract_quasi_rh_at_tau(3.6)
    assert r_quasi["regime"] == "REGIME_II_WINDOW_FREE_QUASI_RH_EXTRACTED"
    assert r_quasi["window_free_quasi_rh_valid"] is True
    assert r_quasi["farkas_certified_in_Q"] is True
    assert r_quasi["farkas_exact_rational_lower_bound"] == "4/5"
    # Verify that ALL deep linear-term automorphic hypotheses gain >= 0.30 powers of slack!
    triage = r_quasi["sorry_hypothesis_triage"]
    for h_slack in (
        "H1_archimedean_size",
        "H2_level_correction",
        "H4_form_count_mass",
        "H5_first_lattice_covolume",
        "H6_lubin_tate_divisibility",
    ):
        assert triage[h_slack]["status"] == "SLACK_DEGRADABLE"
        assert triage[h_slack]["allowable_exponent_loss"] >= 0.299

    r_imp = qrh.extract_quasi_rh_at_tau(2.7)
    assert r_imp["regime"] == "REGIME_III_FARKAS_IMPOSSIBILITY_NO_POWER_SAVING"
    assert r_imp["window_free_quasi_rh_valid"] is False
    assert r_imp["farkas_certified_in_Q"] is True
    assert r_imp["farkas_exact_rational_lower_bound"] == "11/10"
