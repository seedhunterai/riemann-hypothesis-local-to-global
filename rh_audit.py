#!/usr/bin/env python3
"""
rh_audit.py — Local-to-Global Diagnostic & Audit CLI

Built directly on the four use cases separating a one-sided, scale-dependent
RH certificate from an analytic, fixed-object local-to-global derivation:

  1. spectral : Detect double zeros & GUE vs. Poisson statistics that bypass G(x) <= x^(1/2+eps)
  2. curve    : Audit any elliptic curve E/Q (Hasse-Weil, Sato-Tate moments, GL_1 table orthogonality)
  3. table    : Audit Rubinstein's (q, d, N) table vs. exact Dirichlet class-number identity h(-d)
  4. scale    : Quantify (H12) certificate parameter explosion vs. fixed-object effective bounds
  5. report   : Run all four diagnostics in one unified audit

Examples:
    python rh_audit.py report
    python rh_audit.py spectral --num-zeros 60 --inject-double 0
    python rh_audit.py curve --a1 0 --a2 0 --a3 0 --a4 -1 --a6 0 --limit 3000   # CM curve y^2 = x^3 - x
    python rh_audit.py table --q 11 --d 1319 --N 16
    python rh_audit.py scale --x 1e10 --eps 0.01 --json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Any, Dict, List

import mpmath as mp
import numpy as np

import rh_core as rc


# ---------------------------------------------------------------------------
# Subcommand 1: Spectral / Post-RH Blindness Audit (Use Case 1)
# ---------------------------------------------------------------------------

def run_spectral(
    num_zeros: int = 60,
    eps: float = 0.05,
    inject_double: int | None = None,
    poisson: bool = False,
    seed: int = 42,
) -> Dict[str, Any]:
    real_zeros = rc.get_zeta_zeros(num_zeros)
    real_gammas = [z.imag for z in real_zeros]

    if poisson:
        rng = np.random.default_rng(seed)
        mean_gap = (real_gammas[-1] - real_gammas[0]) / (len(real_gammas) - 1)
        gaps = rng.exponential(scale=mean_gap, size=len(real_gammas) - 1)
        gammas = np.concatenate([[real_gammas[0]], real_gammas[0] + np.cumsum(gaps)])
        zeros = [complex(0.5, float(g)) for g in gammas]
        source = "synthetic_poisson_on_critical_line"
    else:
        zeros = real_zeros
        gammas = np.asarray(real_gammas)
        source = "riemann_zeta_zeros"

    envelope_simple = 2.0 + float(mp.euler) - math.log(4.0 * math.pi)
    xs = np.geomspace(1e6, 1e9, 250)

    passes_one_sided_window = True
    max_norm_envelope = 0.0
    for x in xs:
        if inject_double is not None:
            g_val = rc.G_with_double_zero(float(x), zeros, double_idx=inject_double)
        else:
            g_val = rc.G_explicit_formula(float(x), zeros)
        if g_val > float(x) ** (0.5 + eps):
            passes_one_sided_window = False
        norm = abs(g_val - math.log(2.0 * math.pi)) / math.sqrt(float(x))
        if norm > max_norm_envelope:
            max_norm_envelope = norm

    alphas = np.linspace(0.15, 0.35, 25)
    mean_form_factor = float(
        np.mean([rc.pair_correlation_form_factor(gammas, float(a)) for a in alphas])
    )

    multiplicity_flag = "MULTIPLE_ZERO_DETECTED" if max_norm_envelope > envelope_simple * 1.05 else "SIMPLE_CONSISTENT"
    spacing_flag = "GUE_REPULSION" if mean_form_factor < 0.55 else "POISSON_LIKE_UNCORRELATED"
    blindness_triggered = passes_one_sided_window and (
        multiplicity_flag == "MULTIPLE_ZERO_DETECTED" or spacing_flag == "POISSON_LIKE_UNCORRELATED"
    )

    return {
        "source": source,
        "num_zeros": num_zeros,
        "eps": eps,
        "inject_double": inject_double,
        "passes_one_sided_window_bound": passes_one_sided_window,
        "simple_zero_ceiling": round(envelope_simple, 6),
        "observed_max_normalized_G": round(max_norm_envelope, 6),
        "multiplicity_verdict": multiplicity_flag,
        "mean_form_factor_F_0.15_0.35": round(mean_form_factor, 4),
        "spacing_verdict": spacing_flag,
        "one_sided_blindness_triggered": blindness_triggered,
    }


# ---------------------------------------------------------------------------
# Subcommand 2: Elliptic Curve GL_2 vs. GL_1 Audit (Use Case 2)
# ---------------------------------------------------------------------------

def run_curve(
    a1: int = 0,
    a2: int = -1,
    a3: int = 1,
    a4: int = 0,
    a6: int = 0,
    limit: int = 4000,
    q: int = 11,
    d: int = 1319,
) -> Dict[str, Any]:
    primes = [p for p in rc.sieve_primes(limit) if p not in (2, q)]
    traces = []
    hasse_weil_ok = True
    zero_trace_count = 0

    for p in primes:
        ap = rc.elliptic_curve_ap(p, a1=a1, a2=a2, a3=a3, a4=a4, a6=a6)
        if abs(ap) > 2.0 * math.sqrt(p) + 1e-9:
            hasse_weil_ok = False
        if ap == 0:
            zero_trace_count += 1
        traces.append(ap / (2.0 * math.sqrt(p)))

    t = np.asarray(traces, dtype=float)
    m2 = float(np.mean(t ** 2))
    m4 = float(np.mean(t ** 4))
    zero_frac = zero_trace_count / len(primes)

    weights = np.array(
        [(1 + rc.chi_neg(q, p)) * (1 - rc.chi_neg(d, p)) for p in primes],
        dtype=float,
    )
    corr = float(np.corrcoef(t, weights)[0, 1])

    # Classify distribution:
    # - Non-CM elliptic curve (Sato-Tate): M2 ~ 0.25, M4 ~ 0.125, zero_frac ~ 0
    # - CM elliptic curve (over Q): half of primes supersingular (a_p = 0), M2 ~ 0.25, M4 ~ 0.1875, zero_frac ~ 0.5
    if zero_frac > 0.35:
        distribution_verdict = "CM_ELLIPTIC_CURVE (half supersingular, M4 ~ 3/16)"
    elif abs(m2 - 0.25) < 0.05 and abs(m4 - 0.125) < 0.03:
        distribution_verdict = "NON_CM_GL2_SATO_TATE (M2 ~ 1/4, M4 ~ 1/8)"
    else:
        distribution_verdict = "OTHER_OR_SINGULAR"

    return {
        "curve_weierstrass": [a1, a2, a3, a4, a6],
        "primes_tested": len(primes),
        "prime_limit": limit,
        "hasse_weil_holds": hasse_weil_ok,
        "moment_M2": round(m2, 5),
        "moment_M4": round(m4, 5),
        "supersingular_fraction": round(zero_frac, 4),
        "distribution_verdict": distribution_verdict,
        "gl1_table_correlation": round(corr, 5),
        "gl1_table_captures_curve": abs(corr) > 0.2,
    }


# ---------------------------------------------------------------------------
# Subcommand 3: Rubinstein Table vs. Exact Identity Audit (Use Case 3)
# ---------------------------------------------------------------------------

def run_table(
    q: int = 11,
    d: int = 1319,
    N: int = 16,
    a: float = 10000.0,
    w_dist: float = 10.0,
    w_other: float = 0.05,
) -> Dict[str, Any]:
    if N % 4 != 0 or N < 4:
        raise ValueError("N must be a multiple of 4 (N >= 4).")
    n = N // 4
    cols = [2 * i + 1 for i in range(n)]

    lr = rc.L_prime_over_L_at_1(d)
    rows = []
    for c in cols:
        P, M = rc.compute_rubinstein_table(
            q=q, d=d, N=N, a=a, w_dist=w_dist, w_other=w_other, dist_col=c, L_ratio=lr
        )
        rows.append({"col": c, "P": round(P, 3), "M": round(M, 3), "excess": round(P - M, 2)})

    mean_excess = sum(r["excess"] for r in rows) / len(rows)
    max_abs_excess = max(abs(r["excess"]) for r in rows)

    # Exact identity check: h(-d) via reduced forms vs. Dirichlet L(1, chi_{-d})
    h_forms = rc.class_number_imag_quad(d)
    with mp.workdps(30):
        L1 = -sum(
            rc.chi_neg(d, k) * mp.digamma(mp.mpf(k) / d)
            for k in range(1, d)
            if rc.chi_neg(d, k)
        ) / d
        h_analytic = float(L1 * mp.sqrt(d) / mp.pi)

    identity_error = abs(h_analytic - h_forms)
    slack_to_log_h = max_abs_excess / math.log(max(h_forms, 2))

    return {
        "q": q,
        "d": d,
        "N": N,
        "columns": n,
        "L_prime_over_L_1": round(lr, 8),
        "translates": rows,
        "mean_excess": round(mean_excess, 3),
        "max_abs_excess": round(max_abs_excess, 3),
        "class_number_forms": h_forms,
        "class_number_analytic_L1": round(h_analytic, 10),
        "identity_abs_error": float(f"{identity_error:.3e}"),
        "inequality_slack_over_log_h": round(slack_to_log_h, 2),
    }


# ---------------------------------------------------------------------------
# Subcommand 4: Certificate Scale Blowup vs. Fixed-Object Bounds (Use Case 4)
# ---------------------------------------------------------------------------

def run_scale(x: float = 1e10, eps: float = 0.01) -> Dict[str, Any]:
    if x <= 1.0 or eps <= 0.0:
        raise ValueError("Require x > 1 and eps > 0.")

    log10_x = math.log10(x)
    m_exp = 0.5 - eps / 16.0
    k_exp = 4.5 + 10.0 * eps
    n_exp = 2.0 + 8.0 * eps

    err_exponents = {
        "x/m": round(1.0 - m_exp, 6),
        "xm/k": round(1.0 + m_exp - k_exp, 6),
        "m": round(m_exp, 6),
        "k/n^2": round(k_exp - 2.0 * n_exp, 6),
    }
    dominant_exp = max(err_exponents.values())

    log_q = x ** (eps / 5.0)
    log_d = log_q * (math.log(x) ** 2)

    envelope_const = 2.0 + float(mp.euler) - math.log(4.0 * math.pi)
    explicit_bound = envelope_const * math.sqrt(x)

    return {
        "x": x,
        "eps": eps,
        "error_exponents": err_exponents,
        "dominant_error_exponent": dominant_exp,
        "target_exponent_0.5_plus_9eps": round(0.5 + 9.0 * eps, 6),
        "certificate_dimensions_log10": {
            "log10_m_chain_length": round(log10_x * m_exp, 2),
            "log10_n_columns": round(log10_x * n_exp, 2),
            "log10_k_lattice_rank": round(log10_x * k_exp, 2),
            "log10_q_aux_prime": round(log_q / math.log(10.0), 2),
            "log10_d_aux_prime": round(log_d / math.log(10.0), 2),
        },
        "fixed_object_comparison": {
            "elliptic_curve_dim_H1": 2,
            "explicit_formula_G_bound_at_x": round(explicit_bound, 4),
            "computable_explicitly_for_all_x": True,
        },
    }


# ---------------------------------------------------------------------------
# Pretty printers
# ---------------------------------------------------------------------------

def _print_section(title: str, data: Dict[str, Any]) -> None:
    print(f"\n=== {title} ===")
    for k, v in data.items():
        if isinstance(v, list) and v and isinstance(v[0], dict):
            print(f"  {k}:")
            for row in v:
                print("    " + ", ".join(f"{rk}={rv}" for rk, rv in row.items()))
        elif isinstance(v, dict):
            print(f"  {k}:")
            for sk, sv in v.items():
                print(f"    {sk}: {sv}")
        else:
            print(f"  {k}: {v}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rh_audit.py",
        description="Local-to-Global RH Diagnostic & Audit Workbench",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_spec = sub.add_parser("spectral", help="Use Case 1: Double-zero & GUE/Poisson audit")
    p_spec.add_argument("--num-zeros", type=int, default=60)
    p_spec.add_argument("--eps", type=float, default=0.05)
    p_spec.add_argument("--inject-double", type=int, default=None, help="Zero index to make double (e.g. 0)")
    p_spec.add_argument("--poisson", action="store_true", help="Use synthetic Poisson zeros on Re(s)=1/2")
    p_spec.add_argument("--json", action="store_true")

    p_curve = sub.add_parser("curve", help="Use Case 2: Elliptic curve GL_2 vs. GL_1 table audit")
    p_curve.add_argument("--a1", type=int, default=0)
    p_curve.add_argument("--a2", type=int, default=-1)
    p_curve.add_argument("--a3", type=int, default=1)
    p_curve.add_argument("--a4", type=int, default=0)
    p_curve.add_argument("--a6", type=int, default=0)
    p_curve.add_argument("--limit", type=int, default=4000)
    p_curve.add_argument("--q", type=int, default=11)
    p_curve.add_argument("--d", type=int, default=1319)
    p_curve.add_argument("--json", action="store_true")

    p_tab = sub.add_parser("table", help="Use Case 3: Rubinstein table vs. exact h(-d) identity")
    p_tab.add_argument("--q", type=int, default=11)
    p_tab.add_argument("--d", type=int, default=1319)
    p_tab.add_argument("--N", type=int, default=16)
    p_tab.add_argument("--a", type=float, default=10000.0)
    p_tab.add_argument("--w-dist", type=float, default=10.0)
    p_tab.add_argument("--w-other", type=float, default=0.05)
    p_tab.add_argument("--json", action="store_true")

    p_scale = sub.add_parser("scale", help="Use Case 4: Scale-x certificate blowup vs. fixed-object bounds")
    p_scale.add_argument("--x", type=float, default=1e10)
    p_scale.add_argument("--eps", type=float, default=0.01)
    p_scale.add_argument("--json", action="store_true")

    p_rep = sub.add_parser("report", help="Run all 4 diagnostics in one unified audit")
    p_rep.add_argument("--json", action="store_true")

    return parser


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "spectral":
        res = run_spectral(
            num_zeros=args.num_zeros,
            eps=args.eps,
            inject_double=args.inject_double,
            poisson=args.poisson,
        )
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            _print_section("Use Case 1: Spectral / Post-RH Blindness Audit", res)

    elif args.cmd == "curve":
        res = run_curve(
            a1=args.a1, a2=args.a2, a3=args.a3, a4=args.a4, a6=args.a6,
            limit=args.limit, q=args.q, d=args.d,
        )
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            _print_section("Use Case 2: GL_2 Elliptic Curve vs. GL_1 Table Audit", res)

    elif args.cmd == "table":
        res = run_table(
            q=args.q, d=args.d, N=args.N, a=args.a, w_dist=args.w_dist, w_other=args.w_other
        )
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            _print_section("Use Case 3: Inequality vs. Exact Identity Audit", res)

    elif args.cmd == "scale":
        res = run_scale(x=args.x, eps=args.eps)
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            _print_section("Use Case 4: Scale Blowup vs. Fixed-Object Audit", res)

    elif args.cmd == "report":
        full = {
            "use_case_1_spectral_double_zero": run_spectral(num_zeros=40, inject_double=0),
            "use_case_2_curve_E11a1": run_curve(limit=3000),
            "use_case_3_table_d1319": run_table(),
            "use_case_4_scale_1e10": run_scale(x=1e10, eps=0.01),
        }
        if args.json:
            print(json.dumps(full, indent=2))
        else:
            _print_section("Use Case 1: Spectral Audit (Injected Double Zero)", full["use_case_1_spectral_double_zero"])
            _print_section("Use Case 2: GL_2 Curve E_11a1 vs. GL_1 Table", full["use_case_2_curve_E11a1"])
            _print_section("Use Case 3: Rubinstein Table vs. Exact h(-1319)", full["use_case_3_table_d1319"])
            _print_section("Use Case 4: Scale Blowup at x=1e10, eps=0.01", full["use_case_4_scale_1e10"])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
