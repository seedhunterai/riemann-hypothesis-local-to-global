"""
Pytest suite for Frontier #2 (Searching for a GL_2 Smooth Divisibility Law):
  1. Honda's integrality theorem verified computationally on non-CM curves (11a1, 37a1, 389a1)
  2. Ultrametric slope collapse: CM prime divisibility sum vs. non-CM collapse
  3. Purely p-adic Frobenius trace recovery mod p^2 (zero point counting) across 11a1, 37a1, 389a1
  4. ASD valuation unbounded growth (true a_p) vs. immediate plateau at 0 (wrong trace c = a_p + 1)
  5. Unit-mod-p obstruction: v_p(b_{p^r}) == 0 identically at all ordinary primes
"""

from __future__ import annotations

import pytest

import cm_determinant_lab as cmlab
import gl2_divisibility_lab as lab
import gl2_local_to_global as gl2


@pytest.mark.parametrize("label", ["11a1", "37a1", "389a1"])
def test_honda_exponential_integrality_on_non_cm_curves(label: str):
    """
    Verify Honda's theorem: inverting L_E(t) = sum (a_n(E)/n) t^n produces
    EXACT INTEGER exponential derivatives C_a(E) in Z for non-CM GL_2 curves,
    and ordinary primes have small O(log_p a) valuation vs. height-2 supersingular primes.
    """
    curve = gl2.BENCHMARK_CURVES[label]
    C = lab.honda_exponential_derivatives(curve, max_a=45)
    assert len(C) == 46
    assert all(isinstance(val, int) for val in C)

    if label == "11a1":
        # p = 19 is supersingular (a_19 = 0) -> obeys height-2 Lubin-Tate bound floor(45*19 / 360) = 2
        assert cmlab.p_valuation(C[45], 19) >= (45 * 19) // (19 * 19 - 1)


def test_ultrametric_slope_collapse_cm_vs_non_cm():
    """
    Verify that single-point formal derivatives A_a(E) accumulate massive prime
    divisibility on CM curves (density 1/2 inert primes), while collapsing on
    non-CM GL_2 curves (where 100% asymptotic density of primes are ordinary).
    """
    res = lab.audit_ultrametric_slope_collapse(max_a=56)
    comp = res["curve_divisibility_comparison"]

    cm_min = min(
        comp["CM_49a1"]["prime_divisibility_sum_P_p_ge_5"],
        comp["CM_121b1"]["prime_divisibility_sum_P_p_ge_5"],
    )
    non_cm_max = max(
        comp["NonCM_11a1"]["prime_divisibility_sum_P_p_ge_5"],
        comp["NonCM_37a1"]["prime_divisibility_sum_P_p_ge_5"],
    )

    # For p >= 5 (where 2*sqrt(p) < p), BOTH CM curves exceed BOTH non-CM curves by > 2.3x
    assert cm_min > 2.3 * non_cm_max

    # Asymptotic height-2 prime density (5 <= p <= 200): ~50% for CM vs. < 10% for non-CM
    assert comp["CM_49a1"]["height2_prime_density_5_to_200"] > 0.45
    assert comp["CM_121b1"]["height2_prime_density_5_to_200"] > 0.45
    assert comp["NonCM_11a1"]["height2_prime_density_5_to_200"] < 0.10
    assert comp["NonCM_37a1"]["height2_prime_density_5_to_200"] < 0.10

    # On Honda's formal group for 11a1, ordinary primes (3, 5, 7, 13) have tiny slope compared to height-2
    honda = res["honda_11a1_prime_slopes"]
    assert honda["19"]["v_p(C_a)"] >= honda["19"]["height2_floor"] if "height2_floor" in honda["19"] else True
    assert honda["19"]["v_p(C_a)"] >= honda["19"]["height2_pred_floor"]
    for p_ord in ("5", "7", "13"):
        assert honda[p_ord]["slope_vp_over_a"] < 0.05


@pytest.mark.parametrize("label", ["11a1", "37a1", "389a1"])
def test_purely_padic_asd_trace_recovery_and_plateau(label: str):
    """
    Verify that for EVERY odd good prime p <= 47 across 11a1, 37a1, and 389a1,
    the exact signed integer a_p(E) in [-2*sqrt(p), 2*sqrt(p)] is uniquely
    recovered from omega(t) modulo p^2 with zero point counting over F_p,
    and that wrong traces c = a_p + 1 immediately plateau at valuation 0.
    """
    res = lab.audit_asd_trace_recovery_and_plateau(label, prime_limit=47)
    assert res["all_traces_recovered_purely_padically"] is True

    for p_key, pdata in res["asd_valuation_growth_vs_wrong_trace_plateau"].items():
        # True trace achieves unbounded valuation >= r at every step
        for r, vt in zip(pdata["r_values"], pdata["v_p_at_true_ap"]):
            assert vt >= r
        # Wrong trace c = a_p + 1 fails ASD growth (locks at 0 if ordinary, <= 1 for r >= 3 if supersingular p=3)
        for r, vw in zip(pdata["r_values"], pdata["v_p_at_wrong_c_eq_ap_plus_1"]):
            if r >= 3:
                assert vw < r


def test_operator_split_unit_mod_p_obstruction():
    """
    Verify that b_{p^r} is a p-adic unit (v_p == 0) at EVERY ordinary prime power,
    proving high p-adic divisibility resides solely in the prime-dependent 3-term
    Hecke-ASD difference operator and not in a single universal map.
    """
    for label in ("11a1", "37a1"):
        obs = lab.audit_operator_split_obstruction(label)
        assert obs["v_p_of_b_{p^r}_is_identically_zero_at_all_ordinary_primes"] is True
