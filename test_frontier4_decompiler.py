"""
Pytest suite for Frontier #4 (Formal-to-Analytic Proof Decompiler):
  1. Skolem quantifier & constant-choice poset verification (Section 5.2)
  2. Tarjan SCC circularity detection & (H12) scale-uniformity violation catching
  3. Hidden-constraint discovery: Equation (M) alone ((0.5, 1.0, 0.25)) vs.
     full blueprint with (H8)/(H10) codimension constraints ((0.5, 4.5, 2.0))
  4. Exact LP dual shadow prices lambda* = [0.5, 0.0, 0.5, 0.0] pinpointing
     the active bottleneck hypotheses
  5. Exact (H9) Sha-debt slack threshold delta_c = 7/9 and step fragility beta/(beta+1)
"""

from __future__ import annotations

import proof_dag_decompiler as dec


def test_quantifier_poset_well_founded_and_catches_circularity():
    """Verify Section 5.2 acyclicity on the true blueprint and catch injected cycles/leaks."""
    nodes = dec.parse_lean_blueprint(dec.LEAN_RH_BLUEPRINT)
    res = dec.verify_quantifier_poset(nodes)
    assert res["is_well_founded_dag"] is True
    assert res["uniformity_H12_valid"] is True
    assert res["detected_cycles"] == []

    # 1. Inject an illegal scale-uniformity leak: C_unif depends on q_prime (still a DAG, breaks H12!)
    leaked_spec = dec.LEAN_RH_BLUEPRINT.replace(
        "param C_unif    : Real    | depends_on: []",
        "param C_unif    : Real    | depends_on: [q_prime]",
    )
    # Note: x0_thresh depends on C_unif -> Z_window -> x_scale -> Y_aux -> q_prime -> C_unif: forms a true cycle!
    res_cycle = dec.verify_quantifier_poset(dec.parse_lean_blueprint(leaked_spec))
    assert res_cycle["is_well_founded_dag"] is False
    assert res_cycle["uniformity_H12_valid"] is False
    assert len(res_cycle["detected_cycles"]) >= 1

    # 2. Inject a pure non-cyclic uniformity violation (a new constant C_bad downstream of x_scale)
    nodes_leak = dec.parse_lean_blueprint(dec.LEAN_RH_BLUEPRINT)
    for n in nodes_leak:
        if n.name == "C_unif":
            n.depends_on = ["x_scale"]
        if n.name == "x0_thresh":
            n.depends_on = ["eps", "R_centres", "J_root"]  # break cycle so graph stays a DAG
    res_leak = dec.verify_quantifier_poset(nodes_leak)
    assert res_leak["is_well_founded_dag"] is True
    assert res_leak["uniformity_H12_valid"] is False
    assert any("C_unif" in v for v in res_leak["uniformity_violations"])


def test_hidden_constraint_discovery_and_dual_shadow_prices():
    """
    Verify:
      1. Equation (M) alone achieves theta* = 0.5 at minimal exponents (0.5, 1.0, 0.25)
         (3.5 orders of magnitude in x smaller rank k than 4.5!).
      2. Adding (H8)/(H10) codimension constraints recovers Rubinstein's (0.5, 4.5, 2.0).
      3. LP dual shadow prices put exact weight 0.5 on `x/m` and 0.5 on `m` (0 on slack terms).
    """
    analysis = dec.run_hidden_constraint_and_pareto_analysis()

    m_only = analysis["equation_M_alone_minimum_complexity"]
    assert m_only["optimal_abscissa_theta_star"] == 0.5
    assert m_only["minimal_complexity_exponents"] == {
        "alpha_m": 0.5,
        "alpha_k": 1.0,
        "alpha_n": 0.25,
    }

    m_full = analysis["full_blueprint_with_H8_H10_codimension"]
    assert m_full["optimal_abscissa_theta_star"] == 0.5
    assert m_full["minimal_complexity_exponents"] == {
        "alpha_m": 0.5,
        "alpha_k": 4.5,
        "alpha_n": 2.0,
    }

    duals = {t["term"]: t["dual_shadow_price"] for t in m_full["term_breakdown"]}
    assert abs(duals["end_chain (x/m)"] - 0.5) < 1e-5
    assert abs(duals["step_accum (m^1)"] - 0.5) < 1e-5
    assert abs(duals["rank_spread (xm/k * k^0)"]) < 1e-5
    assert abs(duals["table_collision (k/n^2)"]) < 1e-5


def test_sha_debt_buffer_and_step_fragility_pareto():
    """
    Verify:
      1. Because xm/k has 3.5 powers of unused slack at k = x^4.5, full RH (theta* = 0.5)
         survives even if (H9)'s Sha/denominator bound degrades up to k^(7/9) (~ x^3.5),
         and still yields a nontrivial zero-free strip theta* < 1 up to delta = 0.95!
      2. Degrading step accumulation to m^beta matches theta* = beta / (beta + 1) identically.
    """
    analysis = dec.run_hidden_constraint_and_pareto_analysis()

    by_delta = {row["sha_debt_exponent_delta"]: row for row in analysis["H9_sha_debt_slack_pareto_at_k_le_4.5"]}
    assert by_delta[round(7.0 / 9.0, 6)]["still_proves_full_RH"] is True
    assert by_delta[0.85]["still_proves_full_RH"] is False
    assert by_delta[0.85]["proves_nontrivial_zero_free_strip"] is True
    assert by_delta[0.95]["proves_nontrivial_zero_free_strip"] is True

    for row in analysis["step_accumulation_fragility_curve"]:
        assert abs(row["lp_theta_star"] - row["exact_formula_beta_over_beta_plus_1"]) < 1e-5
