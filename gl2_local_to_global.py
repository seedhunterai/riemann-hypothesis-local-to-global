#!/usr/bin/env python3
"""
gl2_local_to_global.py — Frontier GL_2 Local-to-Global Analytic Rank & BSD Solver

Demonstrates the concrete payoff of Use Cases 2, 3, and 4:
Starting ONLY from per-prime local Frobenius point counts a_p = p + 1 - #E(F_p)
of a GL_2 elliptic curve E/Q (conductor N), this tool uses two fixed-object
analytic local-to-global identities that the GL_1 one-sided table cannot formulate:

  1. Exact Central L-Derivative Solver (Cremona/Dokchitser Mellin series):
     Computes L(E,1), L'(E,1), L''(E,1), L'''(E,1) directly from local a_p,
     determining the exact analytic rank r_an(E) and central BSD L-value.
  2. Weil-Mestre Explicit Formula Rank Upper Bound:
     Uses the Odlyzko/Fejer non-negative kernel Phi(i*gamma) >= 0 on local
     Frobenius power sums c(p^k) = alpha_p^k + beta_p^k to certify an explicit
     GRH upper bound on r_an(E).

Examples:
    python gl2_local_to_global.py --benchmark
    python gl2_local_to_global.py --curve 389a1
    python gl2_local_to_global.py --weierstrass 0 1 1 -2 0 --conductor 389 --root-number 1
"""

from __future__ import annotations

import argparse
import functools
import json
import math
from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np

import rh_core as rc


@dataclass(frozen=True)
class EllipticCurveSpec:
    label: str
    a1: int
    a2: int
    a3: int
    a4: int
    a6: int
    conductor: int
    root_number: int  # +1 or -1
    known_rank: int | None = None


BENCHMARK_CURVES: Dict[str, EllipticCurveSpec] = {
    "11a1": EllipticCurveSpec("11a1", 0, -1, 1, 0, 0, conductor=11, root_number=1, known_rank=0),
    "37a1": EllipticCurveSpec("37a1", 0, 0, 1, -1, 0, conductor=37, root_number=-1, known_rank=1),
    "389a1": EllipticCurveSpec("389a1", 0, 1, 1, -2, 0, conductor=389, root_number=1, known_rank=2),
    "5077a1": EllipticCurveSpec("5077a1", 0, 0, 1, -7, 6, conductor=5077, root_number=-1, known_rank=3),
}


# ---------------------------------------------------------------------------
# 1. Local Euler Engine: per-prime a_p -> multiplicative a_n & Frobenius c(p^k)
# ---------------------------------------------------------------------------

def compute_dirichlet_coeffs(curve: EllipticCurveSpec, n_max: int) -> List[int]:
    """
    Compute Dirichlet coefficients a_1, ..., a_{n_max} of L(E, s) = sum a_n n^(-s)
    from local Frobenius traces a_p = p + 1 - #E(F_p).
    """
    primes = rc.sieve_primes(n_max)
    a = [0] * (n_max + 1)
    a[1] = 1

    # Compute prime powers via local Euler factor (1 - a_p T + eps_p p T^2)^(-1)
    for p in primes:
        ap = rc.elliptic_curve_ap(p, curve.a1, curve.a2, curve.a3, curve.a4, curve.a6)
        good = (curve.conductor % p) != 0
        a[p] = ap
        pk = p * p
        prev1 = ap
        prev2 = 1
        while pk <= n_max:
            cur = ap * prev1 - (p * prev2 if good else 0)
            a[pk] = cur
            prev2, prev1 = prev1, cur
            pk *= p

    # Multiplicativity for composite n via smallest prime factor
    spf = list(range(n_max + 1))
    for p in primes:
        if p * p > n_max:
            break
        for m in range(p * p, n_max + 1, p):
            if spf[m] == m:
                spf[m] = p

    for n in range(2, n_max + 1):
        if a[n] != 0:
            continue
        p = spf[n]
        pk = p
        rem = n // p
        while rem % p == 0:
            pk *= p
            rem //= p
        if rem > 1:
            a[n] = a[pk] * a[rem]

    return a


def compute_frobenius_power_sums(
    curve: EllipticCurveSpec, x_max: int
) -> List[Tuple[int, int, float]]:
    """
    Return list of (p^k, p, c(p^k)) for all prime powers p^k <= x_max,
    where c(p^k) = alpha_p^k + beta_p^k at good primes and a_p^k at bad primes.
    """
    primes = rc.sieve_primes(x_max)
    entries: List[Tuple[int, int, float]] = []
    for p in primes:
        ap = rc.elliptic_curve_ap(p, curve.a1, curve.a2, curve.a3, curve.a4, curve.a6)
        good = (curve.conductor % p) != 0
        c_prev2 = 2.0
        c_prev1 = float(ap)
        entries.append((p, p, c_prev1))
        pk = p * p
        while pk <= x_max:
            if good:
                c_cur = ap * c_prev1 - p * c_prev2
            else:
                c_cur = ap * c_prev1
            entries.append((pk, p, c_cur))
            c_prev2, c_prev1 = c_prev1, c_cur
            pk *= p
    return entries


# ---------------------------------------------------------------------------
# 2. Exact Central L-Derivative Solver (Cremona/Dokchitser Mellin series)
# ---------------------------------------------------------------------------

def _completed_L_series(a_coeffs: Sequence[int], N: int, w: int, u: mp.mpf) -> mp.mpf:
    """
    Evaluate the normalized completed L-function:
        L_tilde(1 + u) = (2*pi / sqrt(N)) * Lambda(E, 1 + u)
                       = sum_{n >= 1} (a_n / n) * [x_n^(-u)*Gamma(1+u, x_n) + w*x_n^u*Gamma(1-u, x_n)]
    where x_n = 2*pi*n / sqrt(N).
    Note: L_tilde(1 + u) = (sqrt(N)/(2*pi))^u * Gamma(1 + u) * L(E, 1 + u).
    At the first non-vanishing derivative of order r at u = 0:
        L^(r)(E, 1) = L_tilde^(r)(1)   (since the prefactor equals 1 at u = 0!).
    """
    c = 2 * mp.pi / mp.sqrt(N)
    total = mp.mpf(0)
    for n in range(1, len(a_coeffs)):
        an = a_coeffs[n]
        if an == 0:
            continue
        x = c * n
        if x > 45:
            break
        term1 = (x ** (-u)) * mp.gammainc(1 + u, x)
        term2 = (x ** u) * mp.gammainc(1 - u, x)
        total += (an / mp.mpf(n)) * (term1 + w * term2)
    return total


def solve_analytic_rank_and_L_values(
    curve: EllipticCurveSpec,
    max_order: int = 3,
    zero_tol: float = 1e-6,
) -> Dict[str, Any]:
    """
    Determine the exact analytic rank r_an(E) and leading central derivative
    L^(r)(E, 1) directly from local per-prime Frobenius traces a_p.
    """
    # Need 2*pi*n_max / sqrt(N) >= 35 for ~15-digit convergence
    n_max = max(200, int(math.ceil(35.0 * math.sqrt(curve.conductor) / (2.0 * math.pi))))
    a_coeffs = compute_dirichlet_coeffs(curve, n_max)

    with mp.workdps(30):
        def f(u):
            return _completed_L_series(a_coeffs, curve.conductor, curve.root_number, mp.mpf(u))

        # Compute derivatives of L_tilde(1 + u) at u = 0 using high-precision central stencils
        h = mp.mpf("1e-4")
        f0 = f(0)
        fp1 = f(h)
        fm1 = f(-h)
        fp2 = f(2 * h)
        fm2 = f(-2 * h)

        d0 = f0
        d1 = (fp1 - fm1) / (2 * h)
        d2 = (fp1 - 2 * f0 + fm1) / (h ** 2)
        d3 = (fp2 - 2 * fp1 + 2 * fm1 - fm2) / (2 * (h ** 3))

        derivs = [float(d0), float(d1), float(d2), float(d3)]

    # First non-vanishing derivative gives exact analytic rank r_an and L^(r)(E, 1)
    analytic_rank = None
    leading_value = None
    for r, val in enumerate(derivs[: max_order + 1]):
        if abs(val) > zero_tol:
            analytic_rank = r
            leading_value = val
            break

    return {
        "n_terms_used": n_max,
        "L_tilde_derivatives_at_1": [round(v, 10) for v in derivs],
        "analytic_rank": analytic_rank,
        "leading_L_derivative": round(leading_value, 10) if leading_value is not None else None,
    }


# ---------------------------------------------------------------------------
# 3. Weil-Mestre Explicit Formula GRH Upper Bound on Rank
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=16)
def _archimedean_gamma_term(delta: float) -> float:
    """Archimedean Gamma integral (1/pi) * int Phi(it) * Re(psi(1+it)) dt, cached by delta."""
    with mp.workdps(20):
        d_mp = mp.mpf(delta)

        def integrand(t):
            if t == 0:
                phi = d_mp
                re_psi = -mp.euler
            else:
                u = t * d_mp / 2
                phi = d_mp * (mp.sin(u) / u) ** 2
                re_psi = mp.re(mp.digamma(1 + mp.j * t))
            return (2 / mp.pi) * phi * re_psi

        step = 2 * mp.pi / d_mp
        gamma_term = mp.mpf(0)
        for k in range(25):
            gamma_term += mp.quad(integrand, [k * step, (k + 1) * step])
        return float(gamma_term)


def weil_mestre_rank_bound(
    curve: EllipticCurveSpec,
    delta: float = 7.5,
) -> Dict[str, Any]:
    """
    Compute the Weil-Mestre explicit-formula upper bound on analytic rank r_an(E)
    (under GRH for L(E, s)) using the Odlyzko/Fejer kernel:
        F(t) = (1 - |t|/delta)_+,   Phi(i*gamma) = delta * (sin(gamma*delta/2)/(gamma*delta/2))^2 >= 0.
    """
    x_max = int(math.floor(math.exp(delta)))
    power_sums = compute_frobenius_power_sums(curve, x_max)

    # 1. Conductor term: F(0) * log(N / (2*pi)^2)
    cond_term = math.log(curve.conductor / ((2.0 * math.pi) ** 2))

    # 2. Archimedean Gamma term (cached by delta)
    gamma_term_f = _archimedean_gamma_term(delta)

    # 3. Local prime sum: -2 * sum_{p^k <= e^delta} (c(p^k) * log(p) / p^k) * (1 - log(p^k)/delta)
    prime_term = 0.0
    for pk, p, c_pk in power_sums:
        log_pk = math.log(pk)
        weight = 1.0 - log_pk / delta
        if weight > 0.0:
            prime_term -= 2.0 * (c_pk * math.log(p) / pk) * weight

    rhs = cond_term + gamma_term_f + prime_term
    rank_upper_bound = rhs / delta

    return {
        "delta": delta,
        "prime_cutoff_X": x_max,
        "conductor_term": round(cond_term, 6),
        "archimedean_gamma_term": round(gamma_term_f, 6),
        "local_prime_sum_term": round(prime_term, 6),
        "explicit_formula_RHS": round(rhs, 6),
        "grh_rank_upper_bound": round(rank_upper_bound, 4),
        "certified_integer_rank_ceiling": int(math.floor(rank_upper_bound + 1e-9)),
    }


def analyze_curve(curve: EllipticCurveSpec, delta: float = 7.5) -> Dict[str, Any]:
    exact = solve_analytic_rank_and_L_values(curve)
    weil = weil_mestre_rank_bound(curve, delta=delta)
    return {
        "label": curve.label,
        "weierstrass": [curve.a1, curve.a2, curve.a3, curve.a4, curve.a6],
        "conductor": curve.conductor,
        "root_number": curve.root_number,
        "known_rank": curve.known_rank,
        "exact_local_to_global_solver": exact,
        "weil_mestre_explicit_formula": weil,
        "rank_pinpointed_exactly": (
            exact["analytic_rank"] == curve.known_rank
            and weil["certified_integer_rank_ceiling"] >= exact["analytic_rank"]
        ) if curve.known_rank is not None else None,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gl2_local_to_global.py",
        description="GL_2 Local-to-Global Analytic Rank & Weil-Mestre Explicit Formula Solver",
    )
    parser.add_argument("--benchmark", action="store_true", help="Run on 11a1 (r=0), 37a1 (r=1), 389a1 (r=2), 5077a1 (r=3)")
    parser.add_argument("--curve", choices=sorted(BENCHMARK_CURVES.keys()), default="389a1")
    parser.add_argument("--weierstrass", nargs=5, type=int, metavar=("a1", "a2", "a3", "a4", "a6"))
    parser.add_argument("--conductor", type=int)
    parser.add_argument("--root-number", type=int, choices=[-1, 1])
    parser.add_argument("--delta", type=float, default=7.5, help="Weil-Mestre test-function support delta (X = e^delta)")
    parser.add_argument("--json", action="store_true")

    args = parser.parse_args(argv)

    if args.benchmark:
        results = [analyze_curve(c, delta=args.delta) for c in BENCHMARK_CURVES.values()]
        if args.json:
            print(json.dumps(results, indent=2))
        else:
            print(f"{'Curve':<8} {'N':<6} {'w':<3} {'True r':<7} {'Solved r':<9} {'L^(r)(E,1)':<14} {'Weil Bound':<12} {'Ceiling'}")
            print("-" * 76)
            for r in results:
                ex = r["exact_local_to_global_solver"]
                wm = r["weil_mestre_explicit_formula"]
                print(
                    f"{r['label']:<8} {r['conductor']:<6} {r['root_number']:<+3d} "
                    f"{str(r['known_rank']):<7} {str(ex['analytic_rank']):<9} "
                    f"{ex['leading_L_derivative']:<14.6f} {wm['grh_rank_upper_bound']:<12.4f} "
                    f"<= {wm['certified_integer_rank_ceiling']}"
                )
        return 0

    if args.weierstrass is not None:
        if args.conductor is None or args.root_number is None:
            parser.error("--weierstrass requires --conductor and --root-number")
        a1, a2, a3, a4, a6 = args.weierstrass
        spec = EllipticCurveSpec("custom", a1, a2, a3, a4, a6, args.conductor, args.root_number)
    else:
        spec = BENCHMARK_CURVES[args.curve]

    res = analyze_curve(spec, delta=args.delta)
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
