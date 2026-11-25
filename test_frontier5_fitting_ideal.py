"""
Pytest suite for Frontier #5 (`fitting_ideal_lab.py`):
  1. 100% Padé noise annihilation across ALL 51 shifted scales a0 = 1..201 (up to 324 digits)
     on the Gauss lemniscate via the 3-column Fitting ideal
  2. Multi-rank divided-power Fitting-ideal exactness (r = 1..4, M = r + 2):
     zero super-scale cofactors and scale-independent frozen split-prime factors
  3. General Weierstrass CM curves (d = 7 `49a1`, d = 11 `121b1` at a0 = 80):
     zero split primes, zero super-scale primes, pure Lubin-Tate inert product
  4. GL_2 a_p-free Hecke-Wronskian Fitting ideals across non-CM curves (`11a1`, `37a1`, `389a1`):
     eliminating unknown a_p(E) while achieving v_p >= r across 100% of good primes
"""

from __future__ import annotations

import pytest

import fitting_ideal_lab as lab


def test_lemniscate_fitting_ideal_100pct_pade_annihilation():
    """
    Prove computationally that across ALL 51 scales a0 = 1, 5, ..., 201 (up to 324 digits),
    the 3-column Fitting ideal strips 100% of split primes and 100% of super-scale primes,
    recovering the pure 201-digit height-2 Lubin-Tate inert-prime product!
    """
    res = lab.audit_lemniscate_fitting_ideal(max_a0=201, window_cols=3)
    assert res["scales_tested"] == 51
    assert res["scales_with_any_split_prime"] == 0
    assert res["scales_with_any_super_scale_prime"] == 0
    assert res["100_percent_pade_noise_annihilated"] is True

    top = res["milestones"][-1]  # a0 = 201
    assert top["raw_digits"] == 324
    assert top["fitting_gcd_digits"] == 201
    assert top["pade_noise_digits_stripped"] == 123
    assert top["rel_diff_from_LT_floor_formula"] < 1e-5


def test_multirank_divided_power_fitting_stability():
    """
    Verify that for ranks r = 1, 2, 3, 4 with M = r + 2 columns in the Hurwitz
    divided-power basis, 0 super-scale cofactors survive and split primes freeze
    at a constant independent of scale a0!
    """
    rows = lab.audit_multirank_divided_power_fitting(max_r=4, max_a0=49)
    assert len(rows) == 4
    for row in rows:
        assert row["super_scale_cofactors_found"] == 0
        assert row["split_primes_frozen_independent_of_scale"] is True


def test_general_weierstrass_cm_fitting_ideal_exactness():
    """
    Verify on general Weierstrass CM curves (`49a1` d=7 and `121b1` d=11 at a0 = 80)
    that removing {2, d} from gcd(A_80, A_81, A_82) leaves ZERO split primes and
    ZERO super-scale cofactors (`pure_inert_exactness == True`).
    """
    res = lab.audit_general_cm_fitting_ideal(a0=80)
    for cname in ("CM_d7_49a1", "CM_d11_121b1"):
        assert res[cname]["pure_inert_exactness"] is True
        assert res[cname]["unfactored_super_scale_cofactor"] == 1
        assert res[cname]["split_prime_log_contribution"] == 0.0


@pytest.mark.parametrize("label", ["11a1", "37a1", "389a1"])
def test_gl2_ap_free_hecke_wronskian_divisibility(label: str):
    """
    Verify that the a_p-free 2x2 Hecke-Wronskian Fitting ideal achieves smooth
    valuation v_p >= r across 100% of good primes (ordinary + supersingular)
    on non-CM GL_2 curves without knowing a_p(E) in advance!
    """
    res = lab.audit_gl2_ap_free_wronskian(label)
    assert res["uses_a_p_in_operator"] is False
    assert res["achieves_smooth_divisibility_at_100pct_good_primes"] is True
