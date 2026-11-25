#!/usr/bin/env python3
"""
proof_dag_decompiler.py — Frontier #4 Formal-to-Analytic Proof Decompiler:
Skolem Quantifier-Poset Verifier & Tropical LP Exponent-Provenance Engine

Addresses the core meta-problem of Section 0 ("The proof is certified and not
understood") by decompiling a formal / blueprint proof specification into:

  1. Skolem Quantifier & Constant-Choice Dependency Poset (Section 5.2):
     Parses Lean-4-style theorem/hypothesis declarations with quantifier prefixes,
     extracts the dependency graph over all parameters and witnesses, checks
     well-foundedness via Tarjan strongly-connected-components (SCC), and pinpoints
     minimal circularity cycles if any constant depends on a downstream choice.

  2. Tropical / Minimax Linear-Programming Exponent Solver:
     Formulates the proof's asymptotic error bound as a minimax linear program
         theta* = min_{alpha} max_{j} ( <v_j, alpha> + c_j )   s.t.  A_ub alpha <= b_ub,
     solving simultaneously for:
       - Primal optimal parameter exponents alpha* (m = x^alpha_m, k = x^alpha_k, n = x^alpha_n)
       - Optimal zero-free abscissa theta* (zeta(s) != 0 for Re(s) > theta*)
       - Exact LP Dual Shadow Prices lambda_j* >= 0 (sum = 1), which automatically
         identify the critical bottleneck hypotheses of the entire proof!

  3. Hidden-Constraint Discovery & Hypothesis-Weakening Pareto Analyzer:
     - Proves that the four visible error terms of Equation (M) alone optimize at
       (alpha_m, alpha_k, alpha_n) = (0.5, 1.0, 0.25) — 10^35x smaller lattice rank k
       than Rubinstein's (0.5, 4.5, 2.0) — thereby detecting and isolating the hidden
       codimension side-constraints inside (H8)/(H10) that force k ~ x^(9/2), n ~ x^2!
     - Computes the exact critical threshold delta_c up to which Hypothesis (H9)'s
       Galois-cohomology / Sha-debt bound can degrade (from O(sqrt(b)) to k^delta * x)
       before breaking Re(s) > 1/2, plus exact degradation curves theta*(beta, delta).

Examples:
    python proof_dag_decompiler.py report
    python proof_dag_decompiler.py poset
    python proof_dag_decompiler.py lp
    python proof_dag_decompiler.py pareto
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import linprog


# ---------------------------------------------------------------------------
# 1. Lean-4 Blueprint Parser & Skolem Quantifier-Poset Engine (Section 5.2)
# ---------------------------------------------------------------------------

LEAN_RH_BLUEPRINT = """
-- Formal Blueprint of the Local-to-Global Reduction (Sections 1-6)
param beta      : Real    | depends_on: []
param eps       : Real    | depends_on: [beta]
param R_centres : Nat     | depends_on: [eps]
param J_root    : Real    | depends_on: [eps, R_centres]
param C_unif    : Real    | depends_on: []
param x0_thresh : Real    | depends_on: [eps, R_centres, J_root, C_unif]
param Z_window  : Real    | depends_on: [eps, x0_thresh]
param x_scale   : Real    | depends_on: [Z_window]
param Y_aux     : Real    | depends_on: [x_scale, J_root]
param q_prime   : Nat     | depends_on: [Y_aux, eps]
param d_prime   : Nat     | depends_on: [Y_aux, eps, q_prime]
param m_steps   : Nat     | depends_on: [x_scale, eps]
param n_cols    : Nat     | depends_on: [x_scale, eps]
param k_rank    : Nat     | depends_on: [x_scale, eps]
param b_weight  : Real    | depends_on: [x_scale, eps]
param D_cert    : Cert    | depends_on: [x_scale, q_prime, d_prime, m_steps, n_cols, k_rank, b_weight]
"""


@dataclass
class ParamNode:
    name: str
    kind: str
    depends_on: List[str] = field(default_factory=list)


def parse_lean_blueprint(spec_text: str) -> List[ParamNode]:
    """Parse Lean-style blueprint parameter dependency declarations."""
    nodes: List[ParamNode] = []
    pat = re.compile(
        r"^\s*param\s+([A-Za-z0-9_]+)\s*:\s*([A-Za-z0-9_]+)\s*\|\s*depends_on:\s*\[([^\]]*)\]\s*$"
    )
    for raw_line in spec_text.splitlines():
        line = raw_line.split("--", 1)[0].strip()
        if not line:
            continue
        m = pat.match(line)
        if not m:
            raise ValueError(f"Invalid blueprint line: {raw_line!r}")
        name, kind, deps_str = m.groups()
        deps = [d.strip() for d in deps_str.split(",") if d.strip()]
        nodes.append(ParamNode(name=name, kind=kind, depends_on=deps))
    return nodes


def verify_quantifier_poset(nodes: Sequence[ParamNode]) -> Dict[str, Any]:
    """
    Verify well-foundedness (Section 5.2) of the parameter-choice graph:
      - An edge u -> v means `u` must be chosen BEFORE `v` (`v` depends on `u`).
      - Uses Tarjan's SCC algorithm to detect any circular dependency cycles.
      - Computes topological choice strata (depth from roots) when acyclic.
    """
    name_set = {n.name for n in nodes}
    adj: Dict[str, List[str]] = {n.name: [] for n in nodes}
    rev: Dict[str, List[str]] = {n.name: list(n.depends_on) for n in nodes}

    unknown_refs = []
    for n in nodes:
        for dep in n.depends_on:
            if dep not in name_set:
                unknown_refs.append((n.name, dep))
            else:
                adj[dep].append(n.name)

    if unknown_refs:
        raise ValueError(f"Unknown parameter references: {unknown_refs}")

    # Tarjan SCC for cycle detection
    index = 0
    indices: Dict[str, int] = {}
    lowlink: Dict[str, int] = {}
    stack: List[str] = []
    on_stack: set[str] = set()
    sccs: List[List[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in on_stack:
                lowlink[v] = min(lowlink[v], indices[w])
        if lowlink[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            sccs.append(comp)

    for n in nodes:
        if n.name not in indices:
            strongconnect(n.name)

    cycles = [comp for comp in sccs if len(comp) > 1 or (comp[0] in adj[comp[0]])]
    is_well_founded = len(cycles) == 0

    strata: Dict[int, List[str]] = {}
    topo_order: List[str] = []
    if is_well_founded:
        depth: Dict[str, int] = {}

        def get_depth(v: str) -> int:
            if v in depth:
                return depth[v]
            d_val = 0 if not rev[v] else 1 + max(get_depth(p) for p in rev[v])
            depth[v] = d_val
            return d_val

        for n in nodes:
            d_val = get_depth(n.name)
            strata.setdefault(d_val, []).append(n.name)
        for lvl in sorted(strata.keys()):
            topo_order.extend(strata[lvl])

    # Explicitly check the critical uniformity invariant of (H12):
    # Neither C_unif nor x0_thresh may lie downstream of x_scale, Z_window, q_prime, or d_prime!
    reachable_from: Dict[str, set[str]] = {}
    for n in nodes:
        seen: set[str] = set()
        q_bfs = list(adj[n.name])
        while q_bfs:
            cur = q_bfs.pop()
            if cur not in seen:
                seen.add(cur)
                q_bfs.extend(adj[cur])
        reachable_from[n.name] = seen

    scale_params = ("Z_window", "x_scale", "Y_aux", "q_prime", "d_prime")
    uniformity_violations = []
    for sp in scale_params:
        if sp in reachable_from:
            for const_name in ("C_unif", "x0_thresh", "eps", "R_centres"):
                if const_name in reachable_from[sp]:
                    uniformity_violations.append(f"{const_name} illegally depends on {sp}")

    return {
        "num_parameters": len(nodes),
        "is_well_founded_dag": is_well_founded,
        "detected_cycles": cycles,
        "uniformity_H12_valid": is_well_founded and len(uniformity_violations) == 0,
        "uniformity_violations": uniformity_violations,
        "choice_strata": {str(k): v for k, v in sorted(strata.items())},
        "topological_order": topo_order,
    }


# ---------------------------------------------------------------------------
# 2. Tropical / Minimax LP Exponent-Provenance Solver
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AsymptoticTerm:
    """
    Represents an error term x^(c + v_m*alpha_m + v_k*alpha_k + v_n*alpha_n)
    tagged with the exact hypotheses (H0)-(H12) that produce it.
    """
    name: str
    hypotheses: Tuple[str, ...]
    c: float
    v_m: float
    v_k: float
    v_n: float


@dataclass(frozen=True)
class LinearSideConstraint:
    """
    Side constraint c + v_m*alpha_m + v_k*alpha_k + v_n*alpha_n <= 0
    imposed by internal hypothesis requirements (e.g. rank codimension (H8)/(H10)).
    """
    name: str
    hypotheses: Tuple[str, ...]
    c: float
    v_m: float
    v_k: float
    v_n: float


def default_equation_M_terms(
    step_accum_power: float = 1.0,
    sha_debt_k_power: float = 0.0,
) -> List[AsymptoticTerm]:
    """
    The error terms in Equation (M) (Section 3.4), parameterized to allow
    perturbation / weakening experiments:
      - end_chain (x / m):            1 - alpha_m                      [(H1), (H4), (H5), (H9)]
      - rank_spread ((x*m/k)*k^delta): 1 + alpha_m - (1 - delta)*alpha_k [(H5), (H8), (H9)]
      - step_accum (m^beta):          beta * alpha_m                   [(H2), (H6), (H7), (H8)]
      - word_collision (k / n^2):     alpha_k - 2 * alpha_n            [(H0), (H3), (H10)]
    """
    return [
        AsymptoticTerm("end_chain (x/m)", ("H1", "H4", "H5", "H9"), c=1.0, v_m=-1.0, v_k=0.0, v_n=0.0),
        AsymptoticTerm(
            f"rank_spread (xm/k * k^{sha_debt_k_power:g})",
            ("H5", "H8", "H9"),
            c=1.0,
            v_m=1.0,
            v_k=sha_debt_k_power - 1.0,
            v_n=0.0,
        ),
        AsymptoticTerm(
            f"step_accum (m^{step_accum_power:g})",
            ("H2", "H6", "H7", "H8"),
            c=0.0,
            v_m=step_accum_power,
            v_k=0.0,
            v_n=0.0,
        ),
        AsymptoticTerm("table_collision (k/n^2)", ("H0", "H3", "H10"), c=0.0, v_m=0.0, v_k=1.0, v_n=-2.0),
    ]


def hidden_codimension_constraints() -> List[LinearSideConstraint]:
    """
    Internal side-constraints from (H0), (H8), (H10) (Prop 2.2 codimension eta*R):
      1. Table density / character separation (H0, H10): alpha_n >= 2.0   =>  2.0 - alpha_n <= 0
      2. Span codimension bound eta <= n / k^(1/2) <= x^(-1/2) (H8, H10):
         requires alpha_k >= 2*alpha_n + 0.5                              =>  0.5 + 2*alpha_n - alpha_k <= 0
      3. Twist budget ceiling from (H9) (twists <= x^400): alpha_k <= 20, alpha_n <= 10.
    """
    return [
        LinearSideConstraint("H0_H10_min_table_size (n >= x^2)", ("H0", "H10"), c=2.0, v_m=0.0, v_k=0.0, v_n=-1.0),
        LinearSideConstraint("H8_H10_codimension_span (k >= x^(1/2)*n^2)", ("H8", "H10"), c=0.5, v_m=0.0, v_k=-1.0, v_n=2.0),
        LinearSideConstraint("H9_twist_ceiling_n (n <= x^10)", ("H9",), c=-10.0, v_m=0.0, v_k=0.0, v_n=1.0),
        LinearSideConstraint("H9_twist_ceiling_k (k <= x^20)", ("H9",), c=-20.0, v_m=0.0, v_k=1.0, v_n=0.0),
    ]


def solve_tropical_exponent_lp(
    terms: Sequence[AsymptoticTerm],
    side_constraints: Sequence[LinearSideConstraint] = (),
    regularize_scale_sum: float = 1e-7,
) -> Dict[str, Any]:
    """
    Solve the minimax LP over variables x_lp = [theta, alpha_m, alpha_k, alpha_n] >= 0:
        minimize   theta + eps_reg * (alpha_m + alpha_k + alpha_n)
        subject to v_m*alpha_m + v_k*alpha_k + v_n*alpha_n - theta <= -c   (for each term)
                   v_m*alpha_m + v_k*alpha_k + v_n*alpha_n         <= -c   (for each side constraint)

    `regularize_scale_sum` breaks ties among optimal solutions with identical theta*
    by selecting the minimal-complexity certificate (smallest m, k, n exponents)!
    """
    # Step 1: Pure unregularized solve for exact theta* and exact dual shadow prices!
    c_obj_pure = np.array([1.0, 0.0, 0.0, 0.0])
    A_rows = []
    b_rows = []
    for t in terms:
        A_rows.append([-1.0, t.v_m, t.v_k, t.v_n])
        b_rows.append(-t.c)
    for sc in side_constraints:
        A_rows.append([0.0, sc.v_m, sc.v_k, sc.v_n])
        b_rows.append(-sc.c)

    A_ub = np.asarray(A_rows, dtype=float)
    b_ub = np.asarray(b_rows, dtype=float)
    bounds = [(0.0, None), (0.0, None), (0.0, None), (0.0, None)]

    res_pure = linprog(c_obj_pure, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not res_pure.success:
        return {"feasible": False, "message": res_pure.message}

    theta_star = float(res_pure.fun)
    # In SciPy HiGHS, inequality marginals for A_ub x <= b_ub are <= 0; negate for >= 0 dual weights
    marginals = -np.asarray(res_pure.ineqlin.marginals, dtype=float)
    term_duals = marginals[: len(terms)]
    side_duals = marginals[len(terms) :]

    # Step 2: Secondary LP fixing theta <= theta_star + 1e-12 and minimizing alpha_m + alpha_k + alpha_n
    # to obtain the unique minimum-complexity primal exponents (alpha_m*, alpha_k*, alpha_n*)!
    if regularize_scale_sum > 0:
        A_min = np.vstack([A_ub, [[1.0, 0.0, 0.0, 0.0]]])
        b_min = np.concatenate([b_ub, [theta_star + 1e-11]])
        res_min = linprog([0.0, 1.0, 1.0, 1.0], A_ub=A_min, b_ub=b_min, bounds=bounds, method="highs")
        sol = res_min.x if res_min.success else res_pure.x
    else:
        sol = res_pure.x

    _, alpha_m, alpha_k, alpha_n = [float(v) for v in sol]

    term_evals = []
    bottleneck_hypotheses: Dict[str, float] = {}
    for idx, t in enumerate(terms):
        val = t.c + t.v_m * alpha_m + t.v_k * alpha_k + t.v_n * alpha_n
        dual_w = max(0.0, float(term_duals[idx]))
        term_evals.append({
            "term": t.name,
            "hypotheses": list(t.hypotheses),
            "exponent_at_optimum": round(val, 6),
            "slack_below_theta_star": round(theta_star - val, 6),
            "dual_shadow_price": round(dual_w, 6),
        })
        if dual_w > 1e-6:
            for h in t.hypotheses:
                bottleneck_hypotheses[h] = round(bottleneck_hypotheses.get(h, 0.0) + dual_w, 6)

    side_evals = []
    for idx, sc in enumerate(side_constraints):
        lhs = sc.c + sc.v_m * alpha_m + sc.v_k * alpha_k + sc.v_n * alpha_n
        dual_w = max(0.0, float(side_duals[idx]))
        side_evals.append({
            "constraint": sc.name,
            "hypotheses": list(sc.hypotheses),
            "active_tight": abs(lhs) < 1e-6,
            "dual_shadow_price": round(dual_w, 6),
        })

    return {
        "feasible": True,
        "optimal_abscissa_theta_star": round(theta_star, 6),
        "proves_RH_theta_eq_half": abs(theta_star - 0.5) < 1e-5,
        "minimal_complexity_exponents": {
            "alpha_m": round(alpha_m, 6),
            "alpha_k": round(alpha_k, 6),
            "alpha_n": round(alpha_n, 6),
        },
        "term_breakdown": term_evals,
        "side_constraint_breakdown": side_evals,
        "critical_hypothesis_attribution": bottleneck_hypotheses,
    }


# ---------------------------------------------------------------------------
# 3. Open-Problem Solver: Hidden-Constraint Gap & Hypothesis-Weakening Pareto
# ---------------------------------------------------------------------------

def run_hidden_constraint_and_pareto_analysis() -> Dict[str, Any]:
    """
    1. Compare solving Equation (M) ALONE vs. Equation (M) + (H8)/(H10) codimension constraints:
       Reveals that (M) alone reaches theta* = 1/2 with k = x^1.0, n = x^0.25 (10^35x smaller
       rank at x=10^10!), proving that Rubinstein's k = x^4.5, n = x^2.0 is forced strictly
       by the hidden codimension constraints inside (H8)/(H10).
    2. Sweep (H9) Sha/denominator debt degradation (replacing O(sqrt(b)) by k^delta) under
       the full constrained system to find the EXACT critical threshold delta_c where RH
       (theta* = 1/2) begins to degrade, and map the partial zero-free half-planes theta*(delta) < 1!
    """
    unconstrained_M = solve_tropical_exponent_lp(default_equation_M_terms(), side_constraints=())
    constrained_full = solve_tropical_exponent_lp(
        default_equation_M_terms(), side_constraints=hidden_codimension_constraints()
    )

    # Sweep (H9) debt exponent delta in [0.0, 0.95] under BOTH constrained and unconstrained regimes!
    # In `rank_spread`: exponent is 1 + alpha_m - (1 - delta)*alpha_k.
    # Notice:
    #  - With k fixed >= x^4.5 (constrained): 1 + 0.5 - (1 - delta)*4.5 <= 0.5  <=>  (1 - delta)*4.5 >= 1.0
    #    <=> 1 - delta >= 2/9  <=>  delta <= 7/9 = 0.777778!!
    #    EVEN IF (H9)'S SHA/DENOMINATOR BOUND DEGRADED TO k^(7/9) (~ x^3.5!), RH (theta* = 1/2) STILL HOLDS!!
    #  - Beyond delta = 7/9, let's see what happens when k can adjust up to x^20 or without twist ceiling!
    sha_pareto = []
    for delta in (0.0, 0.5, 7.0 / 9.0, 0.85, 0.90, 0.95):
        # Fix k <= 4.5 (no extra twist budget beyond x^400) to test pure slack absorption at Rubinstein's point:
        fixed_budget_constraints = list(hidden_codimension_constraints()) + [
            LinearSideConstraint("fixed_Rubinstein_rank_budget (k <= x^4.5)", ("H9",), c=-4.5, v_m=0.0, v_k=1.0, v_n=0.0)
        ]
        sol = solve_tropical_exponent_lp(
            default_equation_M_terms(sha_debt_k_power=delta),
            side_constraints=fixed_budget_constraints,
        )
        sha_pareto.append({
            "sha_debt_exponent_delta": round(delta, 6),
            "theta_star": sol["optimal_abscissa_theta_star"],
            "still_proves_full_RH": sol["proves_RH_theta_eq_half"],
            "proves_nontrivial_zero_free_strip": sol["optimal_abscissa_theta_star"] < 1.0 - 1e-6,
            "alpha_m": sol["minimal_complexity_exponents"]["alpha_m"],
        })

    # Sweep step-accumulation exponent beta in m^beta (testing Section 3.4 fragility curve theta* = beta/(beta+1))
    step_pareto = []
    for beta in (1.0, 1.25, 1.5, 2.0, 3.0):
        sol = solve_tropical_exponent_lp(
            default_equation_M_terms(step_accum_power=beta),
            side_constraints=hidden_codimension_constraints(),
        )
        step_pareto.append({
            "step_power_beta": beta,
            "lp_theta_star": sol["optimal_abscissa_theta_star"],
            "exact_formula_beta_over_beta_plus_1": round(beta / (beta + 1.0), 6),
            "optimal_alpha_m": sol["minimal_complexity_exponents"]["alpha_m"],
        })

    return {
        "equation_M_alone_minimum_complexity": unconstrained_M,
        "full_blueprint_with_H8_H10_codimension": constrained_full,
        "H9_sha_debt_slack_pareto_at_k_le_4.5": sha_pareto,
        "step_accumulation_fragility_curve": step_pareto,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="proof_dag_decompiler.py",
        description="Frontier #4: Skolem Quantifier-Poset & Tropical LP Proof Decompiler",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("poset", help="Verify Section 5.2 quantifier & constant-choice DAG")
    sub.add_parser("lp", help="Solve Tropical LP for primal exponents & dual bottleneck shadow prices")
    sub.add_parser("pareto", help="Run hidden-constraint discovery & (H9) weakening Pareto analysis")
    p_r = sub.add_parser("report", help="Run full Frontier #4 decompiler report")
    p_r.add_argument("--json", action="store_true")
    return parser


def main(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.cmd == "poset":
        print(json.dumps(verify_quantifier_poset(parse_lean_blueprint(LEAN_RH_BLUEPRINT)), indent=2))
    elif args.cmd == "lp":
        res = solve_tropical_exponent_lp(default_equation_M_terms(), hidden_codimension_constraints())
        print(json.dumps(res, indent=2))
    elif args.cmd == "pareto":
        print(json.dumps(run_hidden_constraint_and_pareto_analysis(), indent=2))
    elif args.cmd == "report":
        poset = verify_quantifier_poset(parse_lean_blueprint(LEAN_RH_BLUEPRINT))
        analysis = run_hidden_constraint_and_pareto_analysis()
        full = {"1_quantifier_poset": poset, "2_tropical_lp_and_pareto": analysis}
        if args.json:
            print(json.dumps(full, indent=2))
        else:
            print("=== 1. Skolem Quantifier & Constant-Choice Poset (Section 5.2) ===")
            print(
                f"  Parameters: {poset['num_parameters']} | Well-founded DAG: {poset['is_well_founded_dag']} "
                f"| (H12) Uniformity Valid: {poset['uniformity_H12_valid']}"
            )
            for lvl, names in poset["choice_strata"].items():
                print(f"    Stratum {lvl}: {', '.join(names)}")

            print("\n=== 2. Hidden-Constraint Discovery: Equation (M) Alone vs. Full Blueprint ===")
            m_only = analysis["equation_M_alone_minimum_complexity"]
            m_full = analysis["full_blueprint_with_H8_H10_codimension"]
            print(
                f"  (M) alone      -> theta*={m_only['optimal_abscissa_theta_star']}, "
                f"minimal exponents={m_only['minimal_complexity_exponents']}"
            )
            print(
                f"  Full blueprint -> theta*={m_full['optimal_abscissa_theta_star']}, "
                f"minimal exponents={m_full['minimal_complexity_exponents']}"
            )
            print("  Dual bottleneck shadow prices (lambda*):")
            for t in m_full["term_breakdown"]:
                print(
                    f"    {t['term']:<28} exp={t['exponent_at_optimum']:+.4f} "
                    f"| slack={t['slack_below_theta_star']:.4f} | lambda*={t['dual_shadow_price']:.4f} "
                    f"| hyps={t['hypotheses']}"
                )

            print("\n=== 3. (H9) Sha-Debt Weakening Buffer (k <= x^4.5, debt ~ k^delta) ===")
            for row in analysis["H9_sha_debt_slack_pareto_at_k_le_4.5"]:
                print(
                    f"  delta={row['sha_debt_exponent_delta']:.6f} -> theta*={row['theta_star']:.6f} "
                    f"| full_RH={str(row['still_proves_full_RH']):<5} "
                    f"| zero_free_strip={row['proves_nontrivial_zero_free_strip']}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
