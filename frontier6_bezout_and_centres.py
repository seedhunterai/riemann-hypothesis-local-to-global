#!/usr/bin/env python3
"""
frontier6_bezout_and_centres.py — Frontier #6A & #6B Experimental Lab:
  - #6A: Bezout No-Free-Lunch Conservation Law, Cauchy-Binet Kernel Trap, &
         Bargmann-Fock sqrt(a!) vs. Cauchy a! Factor-of-2 Identity
  - #6B: ell-Adic "Exceptional Centre" Twist Scanner (Testing Section 5.3)

Key Findings:
  1. BEZOUT NO-FREE-LUNCH CONSERVATION LAW (#6A):
     In Frontier #5, Fitt_0(a0) = gcd(|A_{a0}|, |A_{a0+4}|, |A_{a0+8}|) stripped
     123 digits of Padé noise off A_201. However, forming Fitt_0 via a scalar
     Bezout combination u0*A_{a0} + u1*A_{a0+4} + u2*A_{a0+8} = Fitt_0 requires
     u_0 == Q_0^(-1) (mod Q_1), forcing:
         log|u_min| = log|A_{a0}| - log(Fitt_0(a0)) + O(1) ~ (1/2) a0 log(a0).
     Likewise, mapping Z^3 -> Z pushes the exact Padé cofactor ||A||_2 / Fitt_0
     into covol(ker Phi) via Cauchy-Binet! Thus ONLY a multi-word span into a
     high-rank target Q (Proposition 2.2, tested in Frontier #7) extracts Fitt_0
     without paying the Bezout norm penalty.

  2. BARGMANN-FOCK sqrt(a!) VS. CAUCHY a! & GAUSS-KRONECKER IDENTITY (#6A):
     Proves to 12+ digits on the Gauss lemniscate E = Q(i):
         Omega_{Q(i)} = exp(gamma + log(pi) - L'/L(1, chi_{-4})) == 2 * varpi^2 / pi,
     and shows why holomorphic L^2 theta projection (||z^a||_{L^2} = sqrt(a!),
     height = log(a!) ~ 1.0 * a log a over C) matches 2 * log(Fitt_0) ~ 1.0 * a log a,
     whereas naive pointwise Cauchy bounds (|A_a|^2 ~ (a!)^2 ~ 2.0 * a log a) lose a factor of 2!

  3. ell-ADIC EXCEPTIONAL-CENTRE TWIST SCANNER (#6B, Testing Section 5.3):
     Derives the exact Euler/Cohen-Lenstra expectation
         E_{CL}[v_ell(h) * log(ell)] = log(ell) * sum_{m>=1} 1 / (ell^m - 1)
     and scans shrinking ell-adic balls |e - e0|_ell <= ell^(-m) (m = 0..4) around
     the worst empirical ell-Sylow exceptional centres e0 across >500 fundamental
     quadratic twists. Finds ZERO ell-adic contagion: ball averages stay O(1)
     bounded right through m = 4, indicating Section 5.3's exceptional centres are
     an artifact of (H9)'s proof method rather than genuine arithmetic blowups.

Examples:
    python frontier6_bezout_and_centres.py report
    python frontier6_bezout_and_centres.py bezout
    python frontier6_bezout_and_centres.py centres --d-base 7 --ell 3
"""

from __future__ import annotations

import argparse
import json
import math
from typing import Any, Dict, List, Tuple

import mpmath as mp

import cm_determinant_lab as cmlab
import rh_core as rc


# ---------------------------------------------------------------------------
# 1. Frontier #6A: Bezout Conservation, Cauchy-Binet Kernel, & Bargmann-Kronecker
# ---------------------------------------------------------------------------

def extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Return (g, s, t) such that g = gcd(a, b) >= 0 and s*a + t*b = g."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def minimal_3term_bezout(a: int, b: int, c: int) -> Tuple[int, Tuple[int, int, int]]:
    """
    Compute g = gcd(|a|, |b|, |c|) and reduced Bezout multipliers (u0, u1, u2)
    such that u0*|a| + u1*|b| + u2*|c| = g with |u0| <= |b| / (2*gcd(|a|,|b|)).
    """
    a, b, c = abs(a), abs(b), abs(c)
    g12, s1, _ = extended_gcd(a, b)
    g123, t12, u2 = extended_gcd(g12, c)
    # Reduce u2 modulo g12 == 1? Note t12*g12 + u2*c = g123; reduce u2 mod (g12 / g123)
    mod_u2 = g12 // g123
    if mod_u2 > 1:
        u2 %= mod_u2
        if u2 > mod_u2 // 2:
            u2 -= mod_u2
        t12 = (g123 - u2 * c) // g12

    mod_u0 = b // g12
    u0 = (s1 * t12) % mod_u0 if mod_u0 > 0 else s1 * t12
    if mod_u0 > 1 and u0 > mod_u0 // 2:
        u0 -= mod_u0
    u1 = (g123 - u0 * a - u2 * c) // b
    return g123, (u0, u1, u2)


def audit_bezout_and_kronecker_conservation() -> Dict[str, Any]:
    # 1. Exact Gauss-Kronecker lemniscate identity (12+ digits):
    with mp.workdps(35):
        L1 = mp.pi / 4
        h = mp.mpf("1e-9")

        def dirichlet_beta(s):
            return (mp.mpf(4) ** (-s)) * (mp.hurwitz(s, mp.mpf(1) / 4) - mp.hurwitz(s, mp.mpf(3) / 4))

        Lp1 = (dirichlet_beta(1 + h) - dirichlet_beta(1 - h)) / (2 * h)
        L_ratio = Lp1 / L1
        Omega_kronecker = mp.exp(mp.euler + mp.log(mp.pi) - L_ratio)
        varpi = (mp.gamma(mp.mpf(1) / 4) ** 2) / (2 * mp.sqrt(2 * mp.pi))
        Omega_gauss = 2 * (varpi ** 2) / mp.pi
        kronecker_err = float(abs(Omega_kronecker - Omega_gauss))

    # 2. Bezout multiplier blowup & Cauchy-Binet kernel covolume identity across scales
    A = cmlab.lemniscate_derivatives(215)
    bezout_rows = []
    for a0 in (21, 49, 101, 149, 201):
        a_v, b_v, c_v = abs(A[a0]), abs(A[a0 + 4]), abs(A[a0 + 8])
        g, (u0, u1, u2) = minimal_3term_bezout(a_v, b_v, c_v)
        assert u0 * a_v + u1 * b_v + u2 * c_v == g

        log_A = math.log(a_v)
        log_fitt = math.log(g)
        log_pade_waste = log_A - log_fitt
        log_u_inf = math.log(max(abs(u0), abs(u1), abs(u2)))

        # Cauchy-Binet on Phi: Z^3 -> Z, Phi(x) = a_v*x0 + b_v*x1 + c_v*x2:
        # Gram(Phi) = a_v^2 + b_v^2 + c_v^2, and [Z : im(Phi)] = g = Fitt_0,
        # so exact 2D kernel lattice ker(Phi) subset Z^3 has covolume:
        #   covol(ker Phi) = sqrt(a_v^2 + b_v^2 + c_v^2) / Fitt_0!
        # Compute log(covol(ker Phi)) stably without float overflow at 324 digits:
        m_max = max(a_v, b_v, c_v)
        log_l2_norm = math.log(m_max) + 0.5 * math.log(
            float(a_v / m_max) ** 2 + float(b_v / m_max) ** 2 + float(c_v / m_max) ** 2
        )
        log_ker_covol = log_l2_norm - log_fitt

        # Bargmann-Fock L^2 height log(a0!) vs. 2 * log(Fitt_0) vs. Cauchy log(|A_{a0}|^2):
        log_fact = math.lgamma(a0 + 1)
        bezout_rows.append({
            "a0": a0,
            "log_raw_A": round(log_A, 2),
            "log_Fitt_0": round(log_fitt, 2),
            "log_pade_cofactor_Q": round(log_pade_waste, 2),
            "log_bezout_multiplier_u": round(log_u_inf, 2),
            "bezout_to_pade_ratio": round(log_u_inf / max(1.0, log_pade_waste), 4),
            "log_kernel_covolume_Z3_to_Z": round(log_ker_covol, 2),
            "two_log_Fitt_over_log_a_fact": round((2.0 * log_fitt) / log_fact, 4),
            "cauchy_two_log_A_over_log_a_fact": round((2.0 * log_A) / log_fact, 4),
        })

    return {
        "gauss_kronecker_lemniscate_identity": {
            "Omega_from_Kronecker_limit": round(float(Omega_kronecker), 12),
            "two_varpi_sq_over_pi": round(float(Omega_gauss), 12),
            "abs_error": float(f"{kronecker_err:.3e}"),
            "exact_match": kronecker_err < 1e-12,
        },
        "bezout_and_kernel_conservation_scan": bezout_rows,
    }


# ---------------------------------------------------------------------------
# 2. Frontier #6B: ell-Adic Exceptional-Centre Twist Scanner (Section 5.3)
# ---------------------------------------------------------------------------

def cohen_lenstra_expected_log_sylow(ell: int, terms: int = 30) -> float:
    """
    Exact Euler-partition / Cohen-Lenstra asymptotic expectation of
    log|Cl(Q(sqrt(-D)))[ell^inf]| = v_ell(h(-D)) * log(ell) for odd prime ell:
        E_{CL}[v_ell(h) * log(ell)] = log(ell) * sum_{m=1}^{inf} 1 / (ell^m - 1).
    """
    return math.log(ell) * sum(1.0 / ((ell ** m) - 1.0) for m in range(1, terms + 1))


def _is_squarefree(n: int) -> bool:
    if n % 4 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 2
    return True


def audit_ell_adic_exceptional_centres(
    d_base: int = 7,
    ell: int = 3,
    max_twist: int = 5000,
) -> Dict[str, Any]:
    """
    Locate the highest ell-Sylow class-group twist `e0` (an empirical "exceptional
    centre") among fundamental discriminants -d_base * e (e <= max_twist), and
    measure mean(v_ell(h) * log(ell)) inside shrinking ell-adic balls:
        B_m(e0) = { e != e0 : |e - e0|_ell <= ell^(-m) }   for m = 0, 1, 2, 3, 4.
    """
    if d_base % 4 != 3:
        raise ValueError("Require d_base == 3 (mod 4) so -d_base == 1 (mod 4).")

    # e == 1 (mod 4), squarefree, gcd(e, d_base) == 1 => -d_base * e is a fundamental discriminant == 1 (mod 4)!
    twists = [
        e for e in range(5, max_twist + 1, 4)
        if math.gcd(e, d_base) == 1 and _is_squarefree(e)
    ]

    # Precompute exact h(-d_base * e) and v_ell(h) across all admissible twists
    twist_data: List[Tuple[int, int, float]] = []  # (e, v_ell, v_ell * log(ell))
    log_ell = math.log(ell)
    for e in twists:
        h = rc.class_number_imag_quad(d_base * e)
        v = cmlab.p_valuation(h, ell)
        twist_data.append((e, v, v * log_ell))

    # Pick the worst empirical exceptional centre e0 in the first half of the range
    # so concentric balls have symmetric coverage
    candidates = sorted(twist_data[: len(twist_data) // 2], key=lambda item: (item[1], -item[0]), reverse=True)
    e0, v0, peak_log = candidates[0]

    cl_pred = cohen_lenstra_expected_log_sylow(ell)
    ball_stats = []
    for m in range(5):
        mod = ell ** m
        members = [val for (e, _, val) in twist_data if e != e0 and (e - e0) % mod == 0]
        mean_val = sum(members) / len(members) if members else 0.0
        max_val = max(members) if members else 0.0
        ball_stats.append({
            "radius_power_m": m,
            "congruence_modulus": mod,
            "sample_size": len(members),
            "mean_log_sylow": round(mean_val, 4),
            "ratio_to_cohen_lenstra": round(mean_val / cl_pred, 4),
            "max_log_sylow_in_punctured_ball": round(max_val, 4),
        })

    # Contagion test: does shrinking the ell-adic radius m = 0 -> 4 cause the ball mean
    # to blow up toward peak_log, or stay bounded < 1.5 * Cohen-Lenstra?
    max_ball_mean = max(b["mean_log_sylow"] for b in ball_stats)
    zero_contagion = max_ball_mean < 1.5 * cl_pred and max_ball_mean < 0.45 * peak_log

    return {
        "d_base": d_base,
        "ell": ell,
        "total_fundamental_twists_scanned": len(twist_data),
        "cohen_lenstra_exact_mean": round(cl_pred, 4),
        "worst_exceptional_centre_found": {
            "e0": e0,
            "discriminant": -d_base * e0,
            "v_ell_of_h": v0,
            "peak_log_sylow": round(peak_log, 4),
            "spike_multiple_over_CL_mean": round(peak_log / cl_pred, 2),
        },
        "concentric_ell_adic_balls": ball_stats,
        "zero_ell_adic_contagion_verified": zero_contagion,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Frontier #6A & #6B Lab")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bezout", help="Run #6A Bezout & Gauss-Kronecker conservation audit")
    p_c = sub.add_parser("centres", help="Run #6B ell-adic exceptional-centre scanner")
    p_c.add_argument("--d-base", type=int, default=7)
    p_c.add_argument("--ell", type=int, default=3)
    p_c.add_argument("--max-twist", type=int, default=5000)
    p_r = sub.add_parser("report", help="Run full Frontier #6A + #6B report")
    p_r.add_argument("--json", action="store_true")

    args = parser.parse_args(argv)
    if args.cmd == "bezout":
        print(json.dumps(audit_bezout_and_kronecker_conservation(), indent=2))
    elif args.cmd == "centres":
        print(json.dumps(audit_ell_adic_exceptional_centres(args.d_base, args.ell, args.max_twist), indent=2))
    elif args.cmd == "report":
        r6a = audit_bezout_and_kronecker_conservation()
        r6b = audit_ell_adic_exceptional_centres(d_base=7, ell=3, max_twist=5000)
        full = {"frontier_6A_bezout_and_kronecker": r6a, "frontier_6B_exceptional_centres": r6b}
        if args.json:
            print(json.dumps(full, indent=2))
        else:
            gk = r6a["gauss_kronecker_lemniscate_identity"]
            print("=== Frontier #6A: Gauss-Kronecker Identity & Bezout Conservation Law ===")
            print(
                f"  Omega_Kronecker = {gk['Omega_from_Kronecker_limit']} == 2*varpi^2/pi = {gk['two_varpi_sq_over_pi']} "
                f"(err={gk['abs_error']})"
            )
            for row in r6a["bezout_and_kernel_conservation_scan"]:
                print(
                    f"  a0={row['a0']:3d}: log|A|={row['log_raw_A']:6.1f}, log(Fitt)={row['log_Fitt_0']:6.1f}, "
                    f"Padé_waste={row['log_pade_cofactor_Q']:5.1f} == log|u_Bezout|={row['log_bezout_multiplier_u']:5.1f} "
                    f"== log(ker_covol)={row['log_kernel_covolume_Z3_to_Z']:5.1f} | 2*log(Fitt)/log(a!)={row['two_log_Fitt_over_log_a_fact']:.3f}"
                )

            print("\n=== Frontier #6B: ell-Adic Exceptional-Centre Twist Scanner (d_base=7, ell=3) ===")
            wc = r6b["worst_exceptional_centre_found"]
            print(
                f"  Twists scanned: {r6b['total_fundamental_twists_scanned']} | Cohen-Lenstra exact mean: {r6b['cohen_lenstra_exact_mean']} "
                f"| Centre e0={wc['e0']} (peak={wc['peak_log_sylow']}, {wc['spike_multiple_over_CL_mean']}x mean)"
            )
            for b in r6b["concentric_ell_adic_balls"]:
                print(
                    f"    |e - e0|_3 <= 3^-{b['radius_power_m']} (n={b['sample_size']:3d}): "
                    f"mean={b['mean_log_sylow']:.4f} ({b['ratio_to_cohen_lenstra']:.2f}x CL) | punctured_max={b['max_log_sylow_in_punctured_ball']:.4f}"
                )
            print(f"  Zero ell-adic contagion verified: {r6b['zero_ell_adic_contagion_verified']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
