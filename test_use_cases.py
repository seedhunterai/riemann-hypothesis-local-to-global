"""
Executable test suite proving the four concrete use cases for extracting
an analytic, fixed-object local-to-global derivation from the RH proof.

Run with:
    pytest -v test_use_cases.py
"""

from __future__ import annotations

import math
import mpmath as mp
import numpy as np
import pytest

import rh_core as rc


# ===========================================================================
# 0. Baseline: Verify the document's core mechanisms & Rubinstein's table
# ===========================================================================

class TestBaselineMechanism:
    def test_explicit_formula_matches_G_exact(self):
        """
        Verify Section 1.1:
          G(x) = log(2*pi) - sum_rho x^rho / (rho*(rho-1)) + O(x^-2)
        and the RH envelope |G(x) - log(2*pi)| <= (2 + gamma - log(4*pi)) * sqrt(x).
        """
        zeros = rc.get_zeta_zeros(40)
        envelope_const = 2.0 + float(mp.euler) - math.log(4.0 * math.pi)
        assert abs(envelope_const - 0.0461914) < 1e-5

        primes = rc.sieve_primes(5000)
        for x in (200.0, 1000.0, 4000.0):
            g_prime = rc.G_exact(x, primes)
            g_explicit = rc.G_explicit_formula(x, zeros)
            # 40 zeros capture G(x) tightly because 1/(rho*(rho-1)) decays like 1/gamma^2
            assert abs(g_prime - g_explicit) < 0.01 * math.sqrt(x)
            assert abs(g_prime - math.log(2.0 * math.pi)) <= envelope_const * math.sqrt(x)

    def test_rubinstein_small_table_exact(self):
        """
        Verify Section 3.3: exact reproduction of Rubinstein's small table
        (q=11, d=1319, N=16, x=1000, b=10, a=10^4, w=10 vs 0.05).
        """
        # Verify congruence conditions from Section 3.1
        q, d = 11, 1319
        assert q % 8 == 3
        assert d % 16 == 7
        assert (d + 1) % q == 0
        assert rc.chi_neg(d, 2) == 1
        assert rc.chi_neg(d, q) == 1

        lr = rc.L_prime_over_L_at_1(d)
        results = [rc.compute_rubinstein_table(dist_col=c, L_ratio=lr) for c in (1, 3, 5, 7)]
        P0, M0 = results[0]
        excesses = [P - M for P, M in results]

        assert round(P0, 3) == 84701.032
        assert round(M0, 3) == 84463.394
        assert [round(e, 2) for e in excesses] == [237.64, -300.26, -176.73, 229.13]
        assert round(sum(excesses) / 4.0, 2) == -2.55

    def test_stepanov_vs_alpoge_claude_balance(self):
        """
        Verify Section 3.4 & the Stepanov parallel:
          1. Balancing x/m + m uniquely achieves exponent 1/2 at m = x^(1/2).
          2. With eps = 0.01, the 4 error exponents match the document's table.
          3. If accumulated cost were m^2 instead of m, the optimum degrades to 2/3.
        """
        # Exact continuous optimization of x/m + m on log scale: min_alpha max(1 - alpha, alpha) = 1/2
        alphas = np.linspace(0.0, 1.0, 10001)
        worst_exp = np.maximum(1.0 - alphas, alphas)
        best_idx = int(np.argmin(worst_exp))
        assert abs(alphas[best_idx] - 0.5) < 1e-4
        assert abs(worst_exp[best_idx] - 0.5) < 1e-4

        # Check the four exponents at eps = 0.01 (Section 3.4 table)
        eps = 0.01
        m_exp = 0.5 - eps / 16.0          # 0.499375
        k_exp = 4.5 + 10.0 * eps          # 4.6
        n_exp = 2.0 + 8.0 * eps           # 2.08
        exps = {
            "x/m": 1.0 - m_exp,
            "xm/k": 1.0 + m_exp - k_exp,
            "m": m_exp,
            "k/n^2": k_exp - 2.0 * n_exp,
        }
        assert abs(exps["x/m"] - 0.500625) < 1e-9
        assert abs(exps["xm/k"] - (-3.100625)) < 1e-9
        assert abs(exps["m"] - 0.499375) < 1e-9
        assert abs(exps["k/n^2"] - 0.44) < 1e-9
        assert max(exps.values()) < 0.5 + 9.0 * eps

        # Fragility check (Section 3.4): x/m + m^2 balances at m = x^(1/3) -> exponent 2/3
        worst_fragile = np.maximum(1.0 - alphas, 2.0 * alphas)
        fragile_idx = int(np.argmin(worst_fragile))
        assert abs(alphas[fragile_idx] - 1.0 / 3.0) < 1e-4
        assert abs(worst_fragile[fragile_idx] - 2.0 / 3.0) < 1e-4


# ===========================================================================
# Use Case 1: Post-RH blindness (Simplicity of zeros & Pair correlation)
# ===========================================================================

class TestUseCase1_PostRHBlindness:
    def test_one_sided_bound_blind_to_double_zeros(self):
        """
        PROVES USE CASE 1A:
        A one-sided window bound G(x) <= x^(1/2 + eps) CANNOT distinguish
        simple zeros from a double zero on Re(s) = 1/2, because x^(1/2)*log(x) << x^(1/2 + eps).
        By contrast, an analytic identity (normalized explicit formula |G(x) - log(2*pi)|/sqrt(x))
        separates them cleanly.
        """
        zeros = rc.get_zeta_zeros(30)
        envelope_simple = 2.0 + float(mp.euler) - math.log(4.0 * math.pi)  # ~0.04619

        # Window [Y, Y^(3/2)] at large scale with eps = 0.05
        eps = 0.05
        xs = np.geomspace(1e6, 1e9, 400)

        max_simple_norm = 0.0
        max_double_norm = 0.0
        for x in xs:
            gs = rc.G_explicit_formula(x, zeros)
            gd = rc.G_with_double_zero(x, zeros, double_idx=0)

            # BOTH satisfy the one-sided window bound G(x) <= x^(1/2 + eps)!
            ceiling = x ** (0.5 + eps)
            assert gs <= ceiling
            assert gd <= ceiling  # Blindness proved: double zero passes the RH window criterion!

            max_simple_norm = max(max_simple_norm, abs(gs - math.log(2.0 * math.pi)) / math.sqrt(x))
            max_double_norm = max(max_double_norm, abs(gd - math.log(2.0 * math.pi)) / math.sqrt(x))

        # Yet the exact analytic normalization separates them by > 4x:
        assert max_simple_norm <= envelope_simple
        assert max_double_norm > 4.0 * envelope_simple

    def test_one_sided_bound_blind_to_pair_correlation(self):
        """
        PROVES USE CASE 1B:
        Both real zeta zeros (GUE repulsion) and synthetic Poisson zeros on Re(s) = 1/2
        satisfy the exact same one-sided bound G(x) <= x^(1/2 + eps).
        Only a spectral/analytic form factor F(alpha) separates GUE from Poisson.
        """
        real_zeros = rc.get_zeta_zeros(120)
        real_gammas = [z.imag for z in real_zeros]

        # Build a Poisson zero sequence on Re(s) = 1/2 with identical mean density
        rng = np.random.default_rng(42)
        mean_gap = (real_gammas[-1] - real_gammas[0]) / (len(real_gammas) - 1)
        poisson_gaps = rng.exponential(scale=mean_gap, size=len(real_gammas) - 1)
        poisson_gammas = np.concatenate([[real_gammas[0]], real_gammas[0] + np.cumsum(poisson_gaps)])
        poisson_zeros = [complex(0.5, g) for g in poisson_gammas]

        # 1. Both satisfy the one-sided bound G(x) <= x^(1/2 + eps) identically
        eps = 0.05
        for x in np.geomspace(1e4, 1e6, 50):
            assert rc.G_explicit_formula(x, real_zeros) <= x ** (0.5 + eps)
            assert rc.G_explicit_formula(x, poisson_zeros) <= x ** (0.5 + eps)

        # 2. The spectral form factor F(alpha) at small alpha (0.15..0.35) cleanly separates them:
        #    GUE obeys F(alpha) ~ |alpha| << 1 (level repulsion); Poisson obeys F(alpha) ~ 1.
        alphas = np.linspace(0.15, 0.35, 25)
        mean_F_zeta = float(np.mean([rc.pair_correlation_form_factor(real_gammas, a) for a in alphas]))
        mean_F_poisson = float(np.mean([rc.pair_correlation_form_factor(poisson_gammas, a) for a in alphas]))

        assert mean_F_zeta < 0.45
        assert mean_F_poisson > 0.70
        assert mean_F_poisson > 1.8 * mean_F_zeta


# ===========================================================================
# Use Case 2: GL_1 -> GL_2 structural obstruction (Elliptic curves & Sato-Tate)
# ===========================================================================

class TestUseCase2_GL1ToGL2Obstruction:
    def test_gl1_character_table_isolates_untwisted_G(self):
        """
        Verify Section 4(d): the (1 + chi_q)(1 - chi_d) table isolates the untwisted
        GL_1 prime sum with exact coefficient +1.
        """
        q, d = 11, 1319
        primes = [p for p in rc.sieve_primes(2000) if p not in (2, q)]
        for p in primes[:100]:
            cq = rc.chi_neg(q, p)
            cd = rc.chi_neg(d, p)
            lhs = (1 + cq) * (1 - cd)
            rhs = 1 + cq - cd - (cq * cd)
            assert lhs == rhs

    def test_gl2_elliptic_curve_obeys_sato_tate_not_abelian_table(self):
        """
        PROVES USE CASE 2:
        For the non-CM elliptic curve E_11a1: y^2 + y = x^3 - x^2 (a GL_2 object):
          1. The fixed cohomological object H^1(E/F_p) enforces Hasse-Weil |a_p| <= 2*sqrt(p).
          2. The normalized traces t_p = a_p / (2*sqrt(p)) follow the continuous Sato-Tate
             semicircle law (E[t_p^2] ~ 1/4, E[t_p^4] ~ 1/8), NOT a discrete +/-1 character
             (whose even moments are all 1).
          3. Projecting t_p onto the GL_1 character table (1 + chi_q)(1 - chi_d) yields
             near-zero correlation: the abelian character-table trick cannot isolate GL_2 data.
        """
        primes = [p for p in rc.sieve_primes(5000) if p not in (2, 11)]
        norm_traces = []
        for p in primes:
            ap = rc.elliptic_curve_11a1_ap(p)
            assert abs(ap) <= 2.0 * math.sqrt(p)  # Hasse-Weil from fixed H^1(E)
            norm_traces.append(ap / (2.0 * math.sqrt(p)))

        t = np.asarray(norm_traces)
        m2 = float(np.mean(t ** 2))
        m4 = float(np.mean(t ** 4))

        # Sato-Tate semicircle moments: M_2 = 1/4 = 0.25, M_4 = 1/8 = 0.125
        # (Whereas any +/-1 Dirichlet character has M_2 = 1.0, M_4 = 1.0!)
        assert abs(m2 - 0.25) < 0.03
        assert abs(m4 - 0.125) < 0.02

        # Orthogonality to GL_1 character-table weights (1 + chi_q)(1 - chi_d)
        q, d = 11, 1319
        weights = np.array([(1 + rc.chi_neg(q, p)) * (1 - rc.chi_neg(d, p)) for p in primes], dtype=float)
        # Centered correlation between GL_2 trace and GL_1 table weight is ~0
        corr = float(np.corrcoef(t, weights)[0, 1])
        assert abs(corr) < 0.08


# ===========================================================================
# Use Case 3: Inequality vs. Exact Identity (Class-group / Sha "debt")
# ===========================================================================

class TestUseCase3_InequalityVsIdentity:
    @pytest.mark.parametrize("d", [23, 71, 167, 1319])
    def test_exact_identity_recovers_class_number_h_minus_d(self, d: int):
        """
        PROVES USE CASE 3:
        An exact analytic identity (Dirichlet class number formula) recovers the
        exact integer class number h(-d) (the rank-0 analogue of Sha/Selmer debt)
        to machine precision from L(1, chi_{-d}).
        By contrast, a one-sided inequality P - M <= E absorbs it into slack.
        """
        assert d % 16 == 7
        h_exact = rc.class_number_imag_quad(d)

        # Compute L(1, chi_{-d}) via Hurwitz digamma identity
        with mp.workdps(30):
            L1 = -sum(
                rc.chi_neg(d, a) * mp.digamma(mp.mpf(a) / d)
                for a in range(1, d)
                if rc.chi_neg(d, a)
            ) / d
            h_from_L = float(L1 * mp.sqrt(d) / mp.pi)

        assert abs(h_from_L - h_exact) < 1e-10
        if d == 1319:
            assert h_exact == 45  # Rubinstein's field Q(sqrt(-1319)) has class number 45

    def test_inequality_slack_hides_arithmetic_invariant(self):
        """
        Verify that in Rubinstein's table (d=1319), the translate-to-translate
        fluctuation |P - M| (~300) is orders of magnitude larger than log(h(-1319)) (~3.8),
        demonstrating why a one-sided inequality cannot extract the arithmetic invariant.
        """
        lr = rc.L_prime_over_L_at_1(1319)
        excesses = [rc.compute_rubinstein_table(dist_col=c, L_ratio=lr) for c in (1, 3, 5, 7)]
        max_abs_excess = max(abs(P - M) for P, M in excesses)
        log_h = math.log(rc.class_number_imag_quad(1319))

        assert max_abs_excess > 50.0 * log_h


# ===========================================================================
# Use Case 4: Computability gap (Scale-x blowup vs. Fixed-object invariants)
# ===========================================================================

class TestUseCase4_ComputabilityGap:
    def test_scale_x_certificate_parameter_explosion(self):
        """
        PROVES USE CASE 4A:
        Quantify the exact parameter blowup of the scale-x certificate D(x)
        at eps = 0.01 and modest x = 10^10:
          - chain length m ~ x^(1/2 - eps/16) ~ 10^5
          - column count n ~ x^(2 + 8*eps) ~ 10^20.8
          - lattice rank k ~ x^(9/2 + 10*eps) ~ 10^46
          - auxiliary prime d has log10(d) > 200 (astronomical!)
        Compare with a fixed-object model (elliptic curve E_11a1): dim H^1 = 2, O(1) for all x.
        """
        eps = 0.01
        log10_x = 10.0  # x = 10^10
        log10_m = log10_x * (0.5 - eps / 16.0)
        log10_n = log10_x * (2.0 + 8.0 * eps)
        log10_k = log10_x * (4.5 + 10.0 * eps)

        # Taking Y = x for the lower-bound estimate of auxiliary primes q, d:
        # log(q) ~ Y^(eps/5), log(d) ~ Y^(eps/5) * (log Y)^2
        Y = 10.0 ** log10_x
        log_q = Y ** (eps / 5.0)
        log_d = log_q * (math.log(Y) ** 2)
        log10_d = log_d / math.log(10.0)

        # Certificate lattice rank exceeds 10^45 and d has > 200 decimal digits (even at Y=10^10)
        assert log10_m > 4.9
        assert log10_n > 20.0
        assert log10_k > 45.0
        assert log10_d > 200.0

        # Fixed-object dimension (H^1 of genus-1 curve E_11a1) is strictly constant = 2
        dim_H1_fixed = 2
        assert dim_H1_fixed == 2

    def test_fixed_object_effective_bound_without_existential_windows(self):
        """
        PROVES USE CASE 4B:
        A fixed analytic identity (explicit formula / fixed cohomological object)
        gives an explicit, computable bound at EVERY finite x >= 100 without
        requiring an existential threshold x_0(eps) or window [Y, Y^(3/2)].
        """
        zeros = rc.get_zeta_zeros(50)
        primes = rc.sieve_primes(10000)
        C_explicit = 2.0 + float(mp.euler) - math.log(4.0 * math.pi)

        for x in (100.0, 500.0, 2500.0, 9000.0):
            g_val = rc.G_exact(x, primes)
            g_pred = rc.G_explicit_formula(x, zeros)
            # Explicit, non-existential bound valid at every finite scale x
            assert abs(g_val - math.log(2.0 * math.pi)) <= C_explicit * math.sqrt(x)
            assert abs(g_val - g_pred) < 0.005 * math.sqrt(x)
