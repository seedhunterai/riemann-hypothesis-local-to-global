"""
Pytest verification suite for Frontiers #6A, #6B, and #7:
  - Frontier #6A (`frontier6_bezout_and_centres.py`):
      * 12-digit Gauss-Kronecker lemniscate limit identity Omega_{Q(i)} == 2*varpi^2/pi
      * Bezout No-Free-Lunch Conservation Law: log|u_min| ~ log|A| - log(Fitt_0)
      * Cauchy-Binet kernel covolume trap: log(covol(ker Phi)) ~ log|A| - log(Fitt_0)
      * Bargmann-Fock sqrt(a!) vs. Cauchy a! factor-of-2 asymptotics
  - Frontier #6B (`frontier6_bezout_and_centres.py`):
      * Exact Euler/Cohen-Lenstra expectation and zero ell-adic contagion across
        shrinking balls |e - e0|_ell <= ell^(-m) around worst Sylow centres (ell=3, 5)
  - Frontier #7 (`frontier7_hnf_slope_engine.py`):
      * Exact modular column Hermite Normal Form (HNF) in rank-deficient regime
        r = D - 1 < D across D = 2, 3, 4: HNF covolume == Fitt_0 identically,
        100% super-scale Padé annihilation, and zero Bezout norm inflation!
"""

from __future__ import annotations

import pytest

import frontier6_bezout_and_centres as f6
import frontier7_hnf_slope_engine as f7


def test_frontier6a_gauss_kronecker_and_bezout_conservation():
    res = f6.audit_bezout_and_kronecker_conservation()

    gk = res["gauss_kronecker_lemniscate_identity"]
    assert gk["exact_match"] is True
    assert gk["abs_error"] < 1e-12

    rows = res["bezout_and_kernel_conservation_scan"]
    top = rows[-1]  # a0 = 201
    # Both Bezout multiplier log-size and Z^3 -> Z kernel covolume grow 1:1 with Padé waste
    # (kernel covolume equals log|A_{a0+8}| - log(Fitt_0), which exceeds log|A_{a0}| - log(Fitt_0) by 8*log(a0/R0))
    assert 0.85 <= top["bezout_to_pade_ratio"] <= 1.15
    assert 1.0 <= top["log_kernel_covolume_Z3_to_Z"] / top["log_pade_cofactor_Q"] <= 1.20

    # Bargmann-Fock 2*log(Fitt_0)/log(a0!) ~ 1.0 vs. Cauchy 2*log|A|/log(a0!) ~ 2.0
    assert 0.95 <= top["two_log_Fitt_over_log_a_fact"] <= 1.25
    assert top["cauchy_two_log_A_over_log_a_fact"] > 1.70


@pytest.mark.parametrize("d_base,ell", [(7, 3), (11, 3), (7, 5)])
def test_frontier6b_zero_ell_adic_exceptional_centre_contagion(d_base: int, ell: int):
    res = f6.audit_ell_adic_exceptional_centres(d_base=d_base, ell=ell, max_twist=4500)
    assert res["worst_exceptional_centre_found"]["spike_multiple_over_CL_mean"] >= 3.5
    assert res["zero_ell_adic_contagion_verified"] is True


@pytest.mark.parametrize("D", [2, 3, 4])
def test_frontier7_rank_deficient_multiword_hnf_exactness(D: int):
    res = f7.simulate_rank_deficient_multiword_hnf(D=D, max_a0=49)
    assert res["each_single_word_is_rank_deficient"] is True
    assert res["hnf_matches_fitting_ideal_at_all_scales"] is True
    assert res["100pct_super_scale_pade_eliminated"] is True

    # At high shifted scales (a0 = 21, 33, 49), single square map suffers > 20% Padé waste
    # while rank-deficient multi-word HNF eliminates 100% of super-scale primes with 0 Bezout penalty!
    by_a0 = {m["a0"]: m for m in res["milestones"]}
    for a0 in (21, 33, 49):
        assert by_a0[a0]["hnf_super_scale_cofactor"] == 1
        assert by_a0[a0]["pade_waste_eliminated_pct"] > 20.0
        assert by_a0[a0]["bezout_multiplier_inflation"] == 0.0
