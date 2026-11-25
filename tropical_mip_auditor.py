#!/usr/bin/env python3
"""
tropical_mip_auditor.py — Stage 2: Disjunctive Tropical MILP & Exact Farkas Auditor

Extends single-cone exponent LP to:
  1. Mixed-Integer Linear Programming (MILP) (`scipy.optimize.milp`) over discrete
     hypothesis choices (e.g., choosing `H9_windowed` vs. `H9_window_free`, or
     selecting the minimum-cost subset of `sorry` hypotheses needed to achieve
     a target zero-free strip Re(s) > theta_target).
  2. Exact Rational Farkas Dual Certificate Generator (`Fraction`-verified):
     Produces a machine-checkable rational dual vector y in Q_{>= 0}^m proving
     algebraically that NO non-negative parameter vector alpha >= 0 can achieve
     theta < theta_lb (certifying impossibility or tightness without floating-point trust).
  3. Bonus Benchmark — Tao-Trudgian-Yang / Guth-Maynard (2024) Zero-Density Polyhedron:
     Solves the reciprocal zero-density polyhedron B(sigma) = 1 / A(sigma) across
     Ingham (1940), Huxley (1972), and Guth-Maynard (2024), certifying both the
     1972 barrier (A* = 12/5 at sigma = 3/4, prime gap 7/12) and the 2024 breakthrough
     (A* = 30/13 at sigma = 7/10, prime gap 17/30) with exact rational dual weights!

Examples:
    python tropical_mip_auditor.py zero-density
    python tropical_mip_auditor.py milp-demo
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

import proof_dag_decompiler as dec


# ---------------------------------------------------------------------------
# 1. Exact Rational Farkas Dual Certificate Extractor & Verifier
# ---------------------------------------------------------------------------

def prove_exact_farkas_lower_bound(
    terms: Sequence[dec.AsymptoticTerm],
    side_constraints: Sequence[dec.LinearSideConstraint] = (),
    max_denom: int = 100_000,
) -> Dict[str, Any]:
    """
    Given primal constraints over x = [theta, alpha_m, alpha_k, alpha_n] >= 0:
        Row_i . x <= b_i   (where Row_i = [-1, v_m, v_k, v_n], b_i = -c_i for terms;
                            Row_i = [ 0, v_m, v_k, v_n], b_i = -c_i for side constraints),
    solve for a non-negative rational dual vector y in Q_{>= 0}^M such that:
        (A^T y)_0 == 1   (wait: with Row_i having -1 on theta, -(A^T y)_0 == 1 => sum_{terms} y_i == 1)
        (A^T y)_j >= 0   for j in {1, 2, 3}  (alpha_m, alpha_k, alpha_n >= 0).
    WAIT — check the sign carefully!
      For every x >= 0 and y >= 0:
          sum_i y_i * (Row_i . x) <= sum_i y_i * b_i.
      Write u = A^T y = [u_theta, u_m, u_k, u_n].
      If u_theta == -1 and u_m >= 0, u_k >= 0, u_n >= 0, then since alpha >= 0:
          -theta <= -theta + u_m*alpha_m + u_k*alpha_k + u_n*alpha_n
                 == u . x == sum_i y_i * (Row_i . x) <= sum_i y_i * b_i.
      Multiplying by -1 proves IDENTICALLY IN EXACT RATIONAL ARITHMETIC:
          theta >= -sum_i y_i * b_i  =  sum_i y_i * c_i !!
    """
    A_rows: List[List[float]] = []
    b_rows: List[float] = []
    labels: List[str] = []

    for t in terms:
        A_rows.append([-1.0, t.v_m, t.v_k, t.v_n])
        b_rows.append(-t.c)
        labels.append(f"term:{t.name}")
    for sc in side_constraints:
        A_rows.append([0.0, sc.v_m, sc.v_k, sc.v_n])
        b_rows.append(-sc.c)
        labels.append(f"constraint:{sc.name}")

    A = np.asarray(A_rows, dtype=float)
    b = np.asarray(b_rows, dtype=float)

    # Dual LP: maximize (-b . y) <=> minimize (b . y) subject to y >= 0,
    #   (A^T y)_0 == -1,  and  -(A^T y)_{1..3} <= 0.
    A_eq = A[:, [0]].T  # shape (1, M)
    b_eq = np.array([-1.0])
    A_ub_dual = -A[:, 1:].T  # shape (3, M): -(A^T y)_{1..3} <= 0  <=>  (A^T y)_{1..3} >= 0
    b_ub_dual = np.zeros(3)

    res = linprog(
        c=b,
        A_ub=A_ub_dual,
        b_ub=b_ub_dual,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=[(0.0, None)] * len(labels),
        method="highs",
    )
    if not res.success:
        return {"certified": False, "message": res.message}

    # Rationalize y and verify 100% in exact `Fraction` arithmetic!
    y_q = [Fraction(float(val)).limit_denominator(max_denom) for val in res.x]

    # Exact rational matrix entries
    A_q: List[List[Fraction]] = []
    b_q: List[Fraction] = []
    for t in terms:
        A_q.append([
            Fraction(-1),
            Fraction(t.v_m).limit_denominator(max_denom),
            Fraction(t.v_k).limit_denominator(max_denom),
            Fraction(t.v_n).limit_denominator(max_denom),
        ])
        b_q.append(Fraction(-t.c).limit_denominator(max_denom))
    for sc in side_constraints:
        A_q.append([
            Fraction(0),
            Fraction(sc.v_m).limit_denominator(max_denom),
            Fraction(sc.v_k).limit_denominator(max_denom),
            Fraction(sc.v_n).limit_denominator(max_denom),
        ])
        b_q.append(Fraction(-sc.c).limit_denominator(max_denom))

    # Verify exact dual feasibility in Q:
    nonneg_ok = all(yi >= 0 for yi in y_q)
    u_q = [sum(y_q[i] * A_q[i][j] for i in range(len(labels))) for j in range(4)]
    exact_valid = nonneg_ok and (u_q[0] == Fraction(-1)) and all(u_q[j] >= 0 for j in (1, 2, 3))

    exact_lb_q = -sum(y_q[i] * b_q[i] for i in range(len(labels)))
    active_weights = {
        labels[i]: str(y_q[i]) for i in range(len(labels)) if y_q[i] > 0
    }

    return {
        "certified": exact_valid,
        "exact_rational_lower_bound": str(exact_lb_q),
        "lower_bound_float": float(exact_lb_q),
        "active_farkas_multipliers": active_weights,
        "dual_column_combination_u": [str(u) for u in u_q],
    }


# ---------------------------------------------------------------------------
# 2. Disjunctive Hypothesis-Selection MILP (`scipy.optimize.milp`)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HypothesisOption:
    """
    A discrete hypothesis variant selectable via binary variable z in {0, 1}:
      - `complexity_cost`: formalization / mathematical weight (e.g. 10.0 for deep Sha/windows, 1.0 for Minkowski)
      - `extra_terms`: asymptotic terms active ONLY when z == 1
      - `extra_constraints`: linear side-constraints active ONLY when z == 1
    """
    option_name: str
    group_name: str
    complexity_cost: float
    uses_windows: bool
    extra_terms: Tuple[dec.AsymptoticTerm, ...] = ()
    extra_constraints: Tuple[dec.LinearSideConstraint, ...] = ()


def solve_hypothesis_selection_milp(
    base_terms: Sequence[dec.AsymptoticTerm],
    base_constraints: Sequence[dec.LinearSideConstraint],
    options: Sequence[HypothesisOption],
    theta_ceiling: float = 0.99,
    big_M: float = 100.0,
) -> Dict[str, Any]:
    """
    Solve a Mixed-Integer Linear Program choosing ONE option per `group_name`
    to MINIMIZE total mathematical/formalization `complexity_cost` (breaking ties
    by minimizing `theta`), subject to achieving a zero-free half-plane:
        theta <= theta_ceiling.
    Variables:
        x = [theta, alpha_m, alpha_k, alpha_n, z_0, z_1, ..., z_{K-1}]
    """
    K = len(options)
    n_vars = 4 + K

    # Objective: minimize sum(cost_l * z_l) + 0.01 * theta
    c_obj = np.zeros(n_vars, dtype=float)
    c_obj[0] = 0.01
    for idx, opt in enumerate(options):
        c_obj[4 + idx] = opt.complexity_cost

    integrality = np.zeros(n_vars, dtype=int)
    integrality[4:] = 1  # binary z_l in {0, 1}

    lb = np.zeros(n_vars, dtype=float)
    ub = np.full(n_vars, np.inf, dtype=float)
    ub[0] = theta_ceiling  # require theta <= theta_ceiling!
    ub[4:] = 1.0

    rows_A: List[List[float]] = []
    rows_lb: List[float] = []
    rows_ub: List[float] = []

    def add_ub_row(vec: List[float], rhs: float) -> None:
        rows_A.append(vec)
        rows_lb.append(-np.inf)
        rows_ub.append(rhs)

    # 1. Base terms & constraints (always active)
    for t in base_terms:
        row = [0.0] * n_vars
        row[0], row[1], row[2], row[3] = -1.0, t.v_m, t.v_k, t.v_n
        add_ub_row(row, -t.c)

    for sc in base_constraints:
        row = [0.0] * n_vars
        row[1], row[2], row[3] = sc.v_m, sc.v_k, sc.v_n
        add_ub_row(row, -sc.c)

    # 2. Exactly-one choice per hypothesis group: sum_{l in G} z_l == 1
    groups: Dict[str, List[int]] = {}
    for idx, opt in enumerate(options):
        groups.setdefault(opt.group_name, []).append(idx)
    for _, idxs in groups.items():
        row = [0.0] * n_vars
        for idx in idxs:
            row[4 + idx] = 1.0
        rows_A.append(row)
        rows_lb.append(1.0)
        rows_ub.append(1.0)

    # 3. Big-M indicator constraints when z_idx == 1:
    #    expr <= -c + big_M * (1 - z_idx)  <=>  expr + big_M * z_idx <= big_M - c
    for idx, opt in enumerate(options):
        for t in opt.extra_terms:
            row = [0.0] * n_vars
            row[0], row[1], row[2], row[3] = -1.0, t.v_m, t.v_k, t.v_n
            row[4 + idx] = big_M
            add_ub_row(row, big_M - t.c)
        for sc in opt.extra_constraints:
            row = [0.0] * n_vars
            row[1], row[2], row[3] = sc.v_m, sc.v_k, sc.v_n
            row[4 + idx] = big_M
            add_ub_row(row, big_M - sc.c)

    constraints = LinearConstraint(np.asarray(rows_A, dtype=float), rows_lb, rows_ub)
    res = milp(c=c_obj, integrality=integrality, bounds=Bounds(lb, ub), constraints=constraints)

    if not res.success:
        return {"feasible": False, "theta_ceiling": theta_ceiling, "message": res.message}

    sol = res.x
    chosen = [options[i] for i in range(K) if sol[4 + i] > 0.5]
    return {
        "feasible": True,
        "theta_ceiling_requested": theta_ceiling,
        "achieved_theta": round(float(sol[0]), 6),
        "exponents": {
            "alpha_m": round(float(sol[1]), 6),
            "alpha_k": round(float(sol[2]), 6),
            "alpha_n": round(float(sol[3]), 6),
        },
        "selected_hypotheses": [c.option_name for c in chosen],
        "requires_windows": any(c.uses_windows for c in chosen),
        "total_complexity_cost": round(sum(c.complexity_cost for c in chosen), 4),
    }


# ---------------------------------------------------------------------------
# 3. Bonus Real-World Target: Guth-Maynard (2024) Zero-Density Polyhedron
# ---------------------------------------------------------------------------

def solve_zero_density_polyhedron(include_guth_maynard: bool = True) -> Dict[str, Any]:
    """
    Solve the Tao-Trudgian-Yang zero-density polyhedron in reciprocal space
    B(sigma) = 1 / A(sigma), where N(sigma, T) << T^(A(sigma)*(1-sigma) + eps):
        B* = min_{sigma in [1/2, 1]} max_{i} B_i(sigma),   A* = 1 / B*,   gap = 1 - B*.
    In reciprocal space, all three historic theorems are EXACT LINEAR FUNCTIONS of sigma:
      - Ingham (1940):       B_Ingham(sigma) = (2 - sigma)/3 = 2/3 - (1/3)*sigma
      - Huxley (1972):       B_Huxley(sigma) = (3*sigma - 1)/3 = sigma - 1/3
      - Guth-Maynard (2024): B_GM(sigma)     = (5*sigma + 3)/15 = (1/3)*sigma + 1/5
    """
    # Variables: x = [B, sigma] >= 0, with 0.5 <= sigma <= 1.0
    # Constraints: B_i(sigma) <= B  <=>  -B + slope_i * sigma <= -intercept_i
    theorems: List[Tuple[str, Fraction, Fraction]] = [
        ("Ingham_1940 (3/(2-sigma))", Fraction(2, 3), Fraction(-1, 3)),
        ("Huxley_1972 (3/(3*sigma-1))", Fraction(-1, 3), Fraction(1, 1)),
    ]
    if include_guth_maynard:
        theorems.append(("Guth_Maynard_2024 (15/(5*sigma+3))", Fraction(1, 5), Fraction(1, 3)))

    A_ub = []
    b_ub = []
    for _, c0, c1 in theorems:
        A_ub.append([-1.0, float(c1)])
        b_ub.append(-float(c0))

    res = linprog(
        c=[1.0, 0.0],
        A_ub=np.asarray(A_ub, dtype=float),
        b_ub=np.asarray(b_ub, dtype=float),
        bounds=[(0.0, None), (0.5, 1.0)],
        method="highs",
    )
    B_float, sigma_float = float(res.x[0]), float(res.x[1])
    duals_float = -np.asarray(res.ineqlin.marginals, dtype=float)

    # Rationalize & verify exact optimality in Q!
    sigma_q = Fraction(sigma_float).limit_denominator(1000)
    B_q = max(c0 + c1 * sigma_q for _, c0, c1 in theorems)
    A_star_q = Fraction(1, 1) / B_q
    gap_q = Fraction(1, 1) - B_q

    duals_q = [Fraction(max(0.0, float(d))).limit_denominator(1000) for d in duals_float]
    # Verify exact Farkas optimality: sum(lambda_i) == 1 and sum(lambda_i * slope_i) == 0!
    dual_sum = sum(duals_q)
    slope_cancel = sum(duals_q[i] * theorems[i][2] for i in range(len(theorems)))
    dual_lb_q = sum(duals_q[i] * theorems[i][1] for i in range(len(theorems)))
    exact_saddle_certified = (
        dual_sum == Fraction(1)
        and slope_cancel == Fraction(0)
        and dual_lb_q == B_q
    )

    return {
        "include_guth_maynard_2024": include_guth_maynard,
        "worst_case_sigma_exact": str(sigma_q),
        "reciprocal_B_star_exact": str(B_q),
        "zero_density_exponent_A_star_exact": str(A_star_q),
        "zero_density_exponent_A_star_float": round(float(A_star_q), 6),
        "prime_gap_exponent_exact": str(gap_q),
        "prime_gap_exponent_float": round(float(gap_q), 6),
        "dual_shadow_prices_exact": {
            theorems[i][0]: str(duals_q[i]) for i in range(len(theorems))
        },
        "exact_rational_saddle_certified": exact_saddle_certified,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 2: Tropical MILP & Zero-Density Polyhedron Auditor")
    parser.add_argument("cmd", choices=["zero-density", "milp-demo"])
    args = parser.parse_args(argv)

    if args.cmd == "zero-density":
        out = {
            "1972_Ingham_Huxley_barrier": solve_zero_density_polyhedron(include_guth_maynard=False),
            "2024_Ingham_Huxley_Guth_Maynard": solve_zero_density_polyhedron(include_guth_maynard=True),
        }
        print(json.dumps(out, indent=2))
    elif args.cmd == "milp-demo":
        opts = [
            HypothesisOption("H9_full_windows_x400", "H9_choice", complexity_cost=10.0, uses_windows=True),
            HypothesisOption(
                "H9_window_free_tau3.6",
                "H9_choice",
                complexity_cost=1.0,
                uses_windows=False,
                extra_constraints=(dec.LinearSideConstraint("tau_k_le_3.6", ("H9_free",), c=-3.6, v_m=0.0, v_k=1.0, v_n=0.0),),
            ),
        ]
        base_terms = dec.default_equation_M_terms() + [
            dec.AsymptoticTerm("table_mod (x^2.5/n)", ("H0",), c=2.5, v_m=0.0, v_k=0.0, v_n=-1.0),
            dec.AsymptoticTerm("codim_span (x*n^2/k)", ("H10",), c=1.0, v_m=0.0, v_k=-1.0, v_n=2.0),
        ]
        out = {
            "min_cost_for_Quasi_RH_theta_0.85": solve_hypothesis_selection_milp(base_terms, (), opts, theta_ceiling=0.85),
            "min_cost_for_Full_RH_theta_0.50": solve_hypothesis_selection_milp(base_terms, (), opts, theta_ceiling=0.50),
        }
        print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
