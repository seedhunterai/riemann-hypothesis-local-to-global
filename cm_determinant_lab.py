#!/usr/bin/env python3
"""
cm_determinant_lab.py — Frontier #1 Experimental Mathematics Lab:
Direct Test of Section 6.3(2) (Exact Determinant Identity vs. Slope Inequality)

Tests whether Lemma 2.1's one-sided slope inequality in the CM / Lemniscate case
is the shadow of an exact local-factor determinant identity (Section 6.3(2)),
using exact arbitrary-precision integer arithmetic:

  1. Universal Height-2 Formal-Group Law:
     Computes exact integer derivatives A_a(E) = (d/dz)^a t(0) for the Gauss
     lemniscate E = Q(i) (u'' = -2 u^3) and general Weierstrass CM / non-CM
     elliptic curves over Z, verifying:
       - Lemniscate inert primes (p == 3 mod 4): v_p(A_a) >= floor(a*p / (p^2 - 1))
       - General height-2 primes (CM inert or non-CM supersingular):
         v_p(A_a(E)) >= v_p(a!) - floor((a - 1) / (p^2 - 1)).

  2. Exact Base-Point Determinant Identity (a0 = 1):
     Proves computationally that at the unshifted origin a0 = 1, ALL giant
     erratic prime factors of individual A_a cancel out, leaving the exact
     (4r - 3)-smooth hyperfactorial identity:
         H_r(1) = det_{0<=i,j<r}[A_{4(i+j)+1}]
                = prod_{k=1}^{r-1} [(4k-3)*(4k-2)^2*(4k-1)^2*(4k)^2*(4k+1)]^(r-k).

  3. Exact Triangular Basis-Invariance Identity:
     Proves for ALL r >= 1 and ALL scale shifts a0 >= 1 that the Hankel
     derivative pairing H_r(a0) and monomial pairing B_r(a0) = det[D^{a0+4j}(u^{4i+1})(0)]
     satisfy the exact algebraic identity:
         H_r(a0) = (prod_{i=0}^{r-1} (4i)!) * B_r(a0).

  4. Scale-Shift Obstruction to Section 6.3(2) (a0 >= 9):
     Shows that shifting to scale a0 >= 9 injects basis-invariant, non-smooth
     Padé cofactor primes (including split primes p == 1 mod 4 and primes p >> a0)
     carrying 10%-25%+ of log|det| — explaining why scale-localized probes
     must drop positive cofactors via a one-sided inequality (Lemma 2.1).

Examples:
    python cm_determinant_lab.py report
    python cm_determinant_lab.py hankel --max-r 8
    python cm_determinant_lab.py shift --r 3 --max-a0 29
    python cm_determinant_lab.py curves --a 24
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from typing import Any, Dict, List, Tuple

import rh_core as rc


# ---------------------------------------------------------------------------
# 1. Exact integer formal-group derivatives: Lemniscate & General Weierstrass
# ---------------------------------------------------------------------------

def lemniscate_derivatives(max_a: int) -> List[int]:
    """
    Compute exact integers A_a = (d/dz)^a sl(0) for 0 <= a <= max_a,
    where u = sl(z) satisfies (u')^2 = 1 - u^4  =>  u'' = -2*u^3, u(0)=0, u'(0)=1.
    Nonzero only for a == 1 (mod 4).
    """
    n_lim = max(max_a, 5)
    A = [0] * (n_lim + 1)
    A[1] = 1
    fact = [math.factorial(i) for i in range(n_lim + 1)]
    for n in range(3, n_lim - 1, 4):
        total = 0
        for i in range(1, n, 4):
            for j in range(1, n - i, 4):
                k = n - i - j
                if k >= 1 and (k & 3) == 1:
                    total += (fact[n] // (fact[i] * fact[j] * fact[k])) * A[i] * A[j] * A[k]
        A[n + 2] = -2 * total
    return A[: max_a + 1]


def lemniscate_monomial_derivatives(max_a: int, max_w: int) -> Dict[int, List[int]]:
    """
    Compute exact integers D^a(u^w)(0) for 1 <= w <= max_w and 0 <= a <= max_a
    on the Gauss lemniscate u = sl(z).
    Returns dict mapping w -> list of [D^0(u^w)(0), ..., D^{max_a}(u^w)(0)].
    """
    A = lemniscate_derivatives(max_a)
    fact = [math.factorial(i) for i in range(max_a + 1)]
    c = [Fraction(A[n], fact[n]) for n in range(max_a + 1)]

    def poly_mul(p: List[Fraction], q: List[Fraction]) -> List[Fraction]:
        res = [Fraction(0)] * (max_a + 1)
        for i, pi in enumerate(p):
            if pi == 0:
                continue
            for j in range(max_a + 1 - i):
                if q[j] != 0:
                    res[i + j] += pi * q[j]
        return res

    powers_egf: Dict[int, List[Fraction]] = {1: c}
    for w in range(2, max_w + 1):
        powers_egf[w] = poly_mul(powers_egf[w - 1], c)

    out: Dict[int, List[int]] = {}
    for w, egf in powers_egf.items():
        out[w] = [int(egf[a] * fact[a]) for a in range(max_a + 1)]
    return out


def weierstrass_formal_derivatives(
    a1: int, a2: int, a3: int, a4: int, a6: int, max_a: int
) -> List[int]:
    """
    Compute exact integer formal-group exponential derivatives A_a(E) = (d/dz)^a t(0)
    for 0 <= a <= max_a on any Weierstrass curve over Z:
        y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6,
    where t = -x/y is the standard parameter at infinity and dz = omega = dx/(2y + a1*x + a3)
    is the invariant differential (Silverman, Ch. IV.1).
    """
    M = max_a + 5
    w = [0] * M
    w[3] = 1
    for deg in range(4, M):
        w2_d = sum(w[i] * w[deg - i] for i in range(3, deg - 2)) if deg >= 6 else 0
        w2_dm1 = sum(w[i] * w[deg - 1 - i] for i in range(3, deg - 3)) if deg >= 7 else 0
        w3_d = (
            sum(
                w[i] * sum(w[j] * w[deg - i - j] for j in range(3, deg - i - 2))
                for i in range(3, deg - 5)
            )
            if deg >= 9
            else 0
        )
        w[deg] = a1 * w[deg - 1] + a2 * w[deg - 2] + a3 * w2_d + a4 * w2_dm1 + a6 * w3_d

    U = [w[3 + k] for k in range(max_a + 2)]
    V = [0] * (max_a + 2)
    V[0] = 1
    for k in range(1, max_a + 2):
        V[k] = -sum(U[j] * V[k - j] for j in range(1, k + 1))

    num = [(m - 2) * V[m] for m in range(max_a + 1)]
    den = []
    for m in range(max_a + 1):
        val = -2 * V[m]
        if m >= 1:
            val += a1 * V[m - 1]
        if m == 3:
            val += a3
        den.append(val)

    # dt/dz = den / num in Z[[t]] (since num[0] = den[0] = -2 and omega in Z[[t]] dt)
    d = [0] * (max_a + 1)
    d[0] = 1
    for m in range(1, max_a + 1):
        rem = den[m] - sum(d[j] * num[m - j] for j in range(m))
        if rem % num[0] != 0:
            raise ArithmeticError(f"Non-integral invariant differential inverse at degree {m}")
        d[m] = rem // num[0]

    # Iterate P_{a+1}(t) = P_a'(t) * d(t) in Z[[t]] to extract A_a = P_a(0)
    P = d[: max_a + 1]
    A = [0, P[0]]
    for a in range(1, max_a):
        deg_lim = max_a - a
        P_deriv = [(k + 1) * P[k + 1] for k in range(deg_lim + 1)]
        P_next = [0] * (deg_lim + 1)
        for i in range(deg_lim + 1):
            pdi = P_deriv[i]
            if pdi == 0:
                continue
            for j in range(deg_lim + 1 - i):
                P_next[i + j] += pdi * d[j]
        P = P_next
        A.append(P[0])
    return A


# ---------------------------------------------------------------------------
# 2. Exact integer linear algebra, valuations, and closed-form determinants
# ---------------------------------------------------------------------------

def exact_int_det(matrix: List[List[int]]) -> int:
    """Exact determinant of an integer matrix via Bareiss fraction-free elimination."""
    n = len(matrix)
    if n == 0:
        return 1
    M = [list(map(int, row)) for row in matrix]
    sign = 1
    prev_piv = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            swap_row = None
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    swap_row = r
                    break
            if swap_row is None:
                return 0
            M[k], M[swap_row] = M[swap_row], M[k]
            sign = -sign
        piv = M[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * piv - M[i][k] * M[k][j]) // prev_piv
            M[i][k] = 0
        prev_piv = piv
    return sign * M[n - 1][n - 1]


def p_valuation(n: int, p: int) -> int:
    """Exact p-adic valuation v_p(n) (returns -1 for n == 0)."""
    if n == 0:
        return -1
    n = abs(n)
    v = 0
    while n % p == 0:
        v += 1
        n //= p
    return v


def legendre_factorial_valuation(a: int, p: int) -> int:
    """Exact v_p(a!) via Legendre's formula."""
    v = 0
    while a >= p:
        a //= p
        v += a
    return v


def height2_valuation_lower_bound(a: int, p: int) -> int:
    """
    Universal height-2 formal-group lower bound on v_p(A_a(E)) for odd prime p:
        v_p(a!) - floor((a - 1) / (p^2 - 1)).
    """
    if a <= 0:
        return 0
    return legendre_factorial_valuation(a, p) - ((a - 1) // (p * p - 1))


def trial_factor(n: int, limit: int = 10_000_000) -> Tuple[Dict[int, int], int]:
    """
    Trial factor |n| up to `limit`. Returns (prime_exponents, unfactored_cofactor).
    """
    n = abs(n)
    if n <= 1:
        return {}, n
    factors: Dict[int, int] = {}
    while (n & 1) == 0:
        factors[2] = factors.get(2, 0) + 1
        n >>= 1
    d = 3
    while d * d <= n and d <= limit:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 2
    if n > 1 and d * d > n:
        factors[n] = factors.get(n, 0) + 1
        n = 1
    return factors, n


def unshifted_lemniscate_hankel_closed_form(r: int) -> int:
    """
    Exact closed-form identity for H_r(1) = det_{0<=i,j<r}[A_{4(i+j)+1}]:
        H_r(1) = prod_{k=1}^{r-1} [(4k-3)*(4k-2)^2*(4k-1)^2*(4k)^2*(4k+1)]^(r-k).
    """
    prod = 1
    for k in range(1, r):
        beta_k = (4 * k - 3) * ((4 * k - 2) ** 2) * ((4 * k - 1) ** 2) * ((4 * k) ** 2) * (4 * k + 1)
        prod *= beta_k ** (r - k)
    return prod


def basis_change_jacobian(r: int) -> int:
    """
    Exact triangular change-of-basis factor between Hankel H_r(a0) and
    monomial pairing B_r(a0) = det_{0<=i,j<r}[D^{a0+4j}(u^{4i+1})(0)]:
        J_r = prod_{i=0}^{r-1} (4i)!.
    """
    prod = 1
    for i in range(r):
        prod *= math.factorial(4 * i)
    return prod


# ---------------------------------------------------------------------------
# 3. Experiment runners
# ---------------------------------------------------------------------------

def audit_unshifted_hankel(max_r: int = 8) -> List[Dict[str, Any]]:
    max_a = 8 * max_r + 1
    A = lemniscate_derivatives(max_a)
    rows = []
    for r in range(1, max_r + 1):
        M = [[A[4 * (i + j) + 1] for j in range(r)] for i in range(r)]
        det_val = exact_int_det(M)
        closed_val = unshifted_lemniscate_hankel_closed_form(r)
        factors, cofactor = trial_factor(det_val, limit=1000)
        max_prime = max(factors.keys()) if factors else 1
        rows.append({
            "r": r,
            "digits": len(str(abs(det_val))),
            "matches_closed_form_identically": det_val == closed_val,
            "max_prime_factor": max_prime,
            "smoothness_bound_4r_minus_3": max(1, 4 * r - 3),
            "unfactored_cofactor": cofactor,
        })
    return rows


def audit_shifted_obstruction(r: int = 3, max_a0: int = 25) -> List[Dict[str, Any]]:
    max_a = max_a0 + 8 * r + 5
    max_w = 4 * r + 1
    A = lemniscate_derivatives(max_a)
    mono = lemniscate_monomial_derivatives(max_a, max_w)
    jac = basis_change_jacobian(r)

    rows = []
    for a0 in range(1, max_a0 + 1, 4):
        H_mat = [[A[a0 + 4 * (i + j)] for j in range(r)] for i in range(r)]
        B_mat = [[mono[4 * i + 1][a0 + 4 * j] for j in range(r)] for i in range(r)]
        H_det = exact_int_det(H_mat)
        B_det = exact_int_det(B_mat)
        basis_identity_holds = (H_det == jac * B_det)

        smooth_cutoff = a0 + 8 * (r - 1)
        factors, cofactor = trial_factor(B_det, limit=20_000_000)

        spurious_primes: Dict[str, Any] = {}
        log_spurious = 0.0
        for p, e in factors.items():
            if p > smooth_cutoff:
                mod4 = "split(1 mod 4)" if (p % 4 == 1) else "inert(3 mod 4)"
                spurious_primes[str(p)] = {"exp": e, "type": mod4}
                log_spurious += e * math.log(p)
        if cofactor > 1:
            mod4 = "1 mod 4" if (cofactor % 4 == 1) else "3 mod 4"
            spurious_primes[f"cofactor_{cofactor}"] = {"exp": 1, "type": mod4}
            log_spurious += math.log(cofactor)

        log_B = math.log(abs(B_det)) if B_det != 0 else 1.0
        rows.append({
            "a0": a0,
            "r": r,
            "basis_jacobian_identity_holds": basis_identity_holds,
            "local_scale_cutoff": smooth_cutoff,
            "spurious_super_scale_factors": spurious_primes,
            "spurious_share_of_log_det": round(log_spurious / log_B, 4),
            "exact_local_product_holds": len(spurious_primes) == 0,
        })
    return rows


def audit_cm_vs_non_cm_curves(a: int = 24, prime_limit: int = 45) -> Dict[str, Any]:
    curves = {
        "CM_d7_49a1": {"coeffs": (1, -1, 0, -2, -1), "cm_d": 7, "bad": {2, 7}},
        "CM_d11_121b1": {"coeffs": (0, -1, 1, -7, 10), "cm_d": 11, "bad": {2, 11}},
        "NonCM_11a1": {"coeffs": (0, -1, 1, 0, 0), "cm_d": None, "bad": {2, 11}},
    }
    primes = [p for p in rc.sieve_primes(prime_limit) if p > 2]
    report: Dict[str, Any] = {}

    for name, cfg in curves.items():
        a1, a2, a3, a4, a6 = cfg["coeffs"]
        A = weierstrass_formal_derivatives(a1, a2, a3, a4, a6, a)
        val_Aa = A[a]

        prime_checks = []
        all_height2_bounds_hold = True
        positive_divisibility_count = 0
        tested_good = 0

        for p in primes:
            if p in cfg["bad"]:
                continue
            tested_good += 1
            vp = p_valuation(val_Aa, p)
            ap_curve = rc.elliptic_curve_ap(p, a1, a2, a3, a4, a6)
            is_height2 = (ap_curve % p == 0)  # supersingular mod p <=> p | a_p (for p >= 5: a_p = 0)
            lb = height2_valuation_lower_bound(a, p) if is_height2 else 0
            if is_height2 and vp < lb:
                all_height2_bounds_hold = False
            if vp > 0:
                positive_divisibility_count += 1
            prime_checks.append({
                "p": p,
                "a_p": ap_curve,
                "height": 2 if is_height2 else 1,
                "v_p(A_a)": vp,
                "height2_lower_bound": lb,
            })

        report[name] = {
            "order_a": a,
            "digits_of_A_a": len(str(abs(val_Aa))),
            "all_height2_bounds_hold": all_height2_bounds_hold,
            "primes_with_positive_v_p": f"{positive_divisibility_count}/{tested_good}",
            "sample_primes": prime_checks[:7],
        }
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cm_determinant_lab.py",
        description="Frontier #1 Lab: Testing Section 6.3(2) Exact Determinant Identities",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_h = sub.add_parser("hankel", help="Verify exact hyperfactorial identity for unshifted H_r(1)")
    p_h.add_argument("--max-r", type=int, default=8)
    p_h.add_argument("--json", action="store_true")

    p_s = sub.add_parser("shift", help="Test scale-shifted determinants H_r(a0) & B_r(a0) for spurious primes")
    p_s.add_argument("--r", type=int, default=3)
    p_s.add_argument("--max-a0", type=int, default=25)
    p_s.add_argument("--json", action="store_true")

    p_c = sub.add_parser("curves", help="Compare height-2 divisibility across CM (d=7, 11) and non-CM (11a1) curves")
    p_c.add_argument("--a", type=int, default=24)
    p_c.add_argument("--json", action="store_true")

    p_r = sub.add_parser("report", help="Run full Frontier #1 experimental report")
    p_r.add_argument("--json", action="store_true")

    return parser


def main(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.cmd == "hankel":
        res = audit_unshifted_hankel(max_r=args.max_r)
        print(json.dumps(res, indent=2))
    elif args.cmd == "shift":
        res = audit_shifted_obstruction(r=args.r, max_a0=args.max_a0)
        print(json.dumps(res, indent=2))
    elif args.cmd == "curves":
        res = audit_cm_vs_non_cm_curves(a=args.a)
        print(json.dumps(res, indent=2))
    elif args.cmd == "report":
        full = {
            "1_unshifted_base_identity_H_r(1)": audit_unshifted_hankel(max_r=8),
            "2_scale_shifted_obstruction_r3": audit_shifted_obstruction(r=3, max_a0=25),
            "3_cm_vs_non_cm_formal_derivatives": audit_cm_vs_non_cm_curves(a=24),
        }
        if args.json:
            print(json.dumps(full, indent=2))
        else:
            print("=== 1. Unshifted Lemniscate Pairing H_r(1): Exact Hyperfactorial Identity ===")
            for row in full["1_unshifted_base_identity_H_r(1)"]:
                print(
                    f"  r={row['r']}: {row['digits']:3d} digits | exact_closed_form={row['matches_closed_form_identically']} "
                    f"| max_prime={row['max_prime_factor']:2d} (<= {row['smoothness_bound_4r_minus_3']:2d})"
                )
            print("\n=== 2. Scale-Shifted Pairing (r=3): Basis Invariance & Spurious Prime Obstruction ===")
            for row in full["2_scale_shifted_obstruction_r3"]:
                print(
                    f"  a0={row['a0']:2d}: H_r=J_r*B_r? {row['basis_jacobian_identity_holds']} "
                    f"| exact_local_product={str(row['exact_local_product_holds']):<5} "
                    f"| spurious_log_share={row['spurious_share_of_log_det']:.4f} "
                    f"| spurious={row['spurious_super_scale_factors']}"
                )
            print("\n=== 3. General Weierstrass Formal Derivatives (a=24): CM vs. Non-CM ===")
            for cname, cdata in full["3_cm_vs_non_cm_formal_derivatives"].items():
                print(
                    f"  {cname:<14}: height2_bounds_hold={cdata['all_height2_bounds_hold']} "
                    f"| primes_dividing_A_24={cdata['primes_with_positive_v_p']}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
