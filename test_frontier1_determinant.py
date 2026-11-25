"""
Pytest suite for Frontier #1 (Section 6.3(2) Direct Experimental Test):
  1. Section 3.2(a) lemniscate inert-prime floor law v_p(A_a) >= floor(a*p/(p^2-1))
  2. Universal Weierstrass height-2 formal-group law across CM (d=7, 11) and non-CM (11a1)
  3. Exact unshifted hyperfactorial determinant identity H_r(1) up to 214 digits (r=8)
  4. Exact basis-change Jacobian identity H_r(a0) = (prod_{i=0}^{r-1} (4i)!) * B_r(a0)
  5. Scale-shift obstruction to Section 6.3(2): emergence of basis-invariant spurious
     split and super-scale primes for a0 >= 9
"""

from __future__ import annotations

import cm_determinant_lab as lab
import rh_core as rc


def test_lemniscate_section_3_2_a_inert_divisibility():
    """Verify Section 3.2(a): v_p(A_a) >= floor(a*p/(p^2-1)) for all inert p == 3 mod 4."""
    max_a = 65
    A = lab.lemniscate_derivatives(max_a)
    inert_primes = [p for p in rc.sieve_primes(70) if p % 4 == 3]

    for a in range(1, max_a + 1, 4):
        assert A[a] != 0
        # Exact 2-adic valuation v_2(A_{4k+1}) = (a - 1) // 2
        assert lab.p_valuation(A[a], 2) == (a - 1) // 2
        for p in inert_primes:
            vp = lab.p_valuation(A[a], p)
            floor_bound = (a * p) // (p * p - 1)
            assert vp >= floor_bound


def test_unshifted_hankel_exact_hyperfactorial_identity():
    """
    Verify that at the base point a0 = 1, Section 6.3(2) holds:
    all giant erratic primes in individual A_a cancel out, matching the exact
    (4r - 3)-smooth hyperfactorial product up to r = 8 (208 digits).
    """
    rows = lab.audit_unshifted_hankel(max_r=8)
    assert rows[-1]["digits"] == 208
    for row in rows:
        assert row["matches_closed_form_identically"] is True
        assert row["unfactored_cofactor"] == 1
        assert row["max_prime_factor"] <= row["smoothness_bound_4r_minus_3"]


def test_basis_change_jacobian_and_shifted_obstruction():
    """
    Verify:
      1. Exact algebraic identity H_r(a0) = (prod_{i=0}^{r-1} (4i)!) * B_r(a0) for all a0.
      2. Unshifted branches a0 in {1, 5} have zero super-scale primes.
      3. Scale-shifted branches a0 >= 9 inject basis-invariant spurious super-scale primes
         (including split primes == 1 mod 4) carrying > 10% of log|det|.
    """
    rows = lab.audit_shifted_obstruction(r=3, max_a0=25)
    by_a0 = {row["a0"]: row for row in rows}

    for row in rows:
        assert row["basis_jacobian_identity_holds"] is True

    # Base Stieltjes branches (a0 = 1, 5) obey exact smooth product
    assert by_a0[1]["exact_local_product_holds"] is True
    assert by_a0[5]["exact_local_product_holds"] is True

    # Shifted branches (a0 >= 9) shatter exact local-factor smoothness
    for a0 in (9, 13, 17, 21, 25):
        assert by_a0[a0]["exact_local_product_holds"] is False
        assert by_a0[a0]["spurious_share_of_log_det"] > 0.10

    # At r=2, a0=17 specifically produces the 6-digit SPLIT prime 131893 (131893 % 4 == 1)
    rows_r2 = {row["a0"]: row for row in lab.audit_shifted_obstruction(r=2, max_a0=21)}
    assert "131893" in rows_r2[17]["spurious_super_scale_factors"]
    assert rows_r2[17]["spurious_super_scale_factors"]["131893"]["type"] == "split(1 mod 4)"


def test_general_weierstrass_height2_law_cm_and_non_cm():
    """
    Verify universal height-2 valuation bound across general Weierstrass curves
    (CM d=7, CM d=11, and non-CM 11a1 at supersingular primes), and confirm CM
    curves have markedly higher density of primes dividing A_a than non-CM 11a1.
    """
    res = lab.audit_cm_vs_non_cm_curves(a=40, prime_limit=75)
    for cname, cdata in res.items():
        assert cdata["all_height2_bounds_hold"] is True

    cm_d11_hits = int(res["CM_d11_121b1"]["primes_with_positive_v_p"].split("/")[0])
    cm_d7_hits = int(res["CM_d7_49a1"]["primes_with_positive_v_p"].split("/")[0])
    non_cm_hits = int(res["NonCM_11a1"]["primes_with_positive_v_p"].split("/")[0])
    assert cm_d11_hits >= 2 * non_cm_hits
    assert cm_d7_hits >= 2 * non_cm_hits
