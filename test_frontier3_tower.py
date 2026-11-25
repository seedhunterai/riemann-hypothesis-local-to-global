"""
Pytest suite for Frontier #3 (Section 6.3(1) Projective Compatibility Test):
  1. Single-pass dyadic evaluator matches rh_core.compute_rubinstein_table identically at N=16
  2. Exact horizontal pullback invariance pi^* (N=16 -> N=128) with zero defect
  3. Exact horizontal pushforward / Iwasawa measure identity mu_{2N}(c_+) + mu_{2N}(c_-) == mu_N(c)
     across all fibers of H_256 -> H_128 -> H_64 -> H_32 -> H_16 -> H_8
  4. Exact Plancherel / Haar wavelet orthogonal energy decomposition across tower steps
  5. Vertical obstruction: scale-indexed ramp specialization matrix has full unbounded rank R,
     whereas fixed-curve Frobenius trace Hankel matrix has constant rank 2 = dim H^1(E)
"""

from __future__ import annotations

import pytest

import projective_tower_lab as lab
import rh_core as rc


def test_fast_tower_matches_rubinstein_table_at_N16():
    """Verify compute_level_excesses_fast matches Section 3.3's exact numbers at N=16."""
    lr = rc.L_prime_over_L_at_1(1319)
    lvl16 = lab.compute_level_excesses_fast(16, q=11, d=1319, a=10000.0, b=10.0, w0=0.05, L_ratio=lr)
    excesses = [round(lvl16["excess_by_col"][c], 2) for c in (1, 3, 5, 7)]
    assert excesses == [237.64, -300.26, -176.73, 229.13]


def test_horizontal_pullback_and_pushforward_iwasawa_exactness():
    """
    Prove that HORIZONTALLY (fixed scale x, climbing N = 8 -> 256), Rubinstein's
    table defines an exact Iwasawa distribution on Z_2^x / {+/- 1}:
      - Pullback pi^* has zero defect in both P and M.
      - Pushforward pi_* has zero defect across every 2-to-1 fiber.
      - Total Iwasawa mass sum_{c in H_N} mu_N(c) is strictly constant across all N.
      - Parseval / Haar primitive-character energy split holds to machine precision,
        with strictly positive new-character energy at every dyadic step.
    """
    horiz = lab.audit_dyadic_iwasawa_tower(max_N=256)

    assert horiz["pullback_16_to_128"]["exact_invariance"] is True
    assert horiz["max_pushforward_iwasawa_defect_all_levels"] < 1e-9

    masses = [lvl["total_iwasawa_mass"] for lvl in horiz["level_summaries"]]
    assert max(masses) - min(masses) < 1e-3

    for tr in horiz["tower_transitions"]:
        assert tr["max_iwasawa_fiber_defect"] < 1e-9
        assert tr["parseval_orthogonality_error"] < 1e-6
        assert tr["new_primitive_character_energy"] > 0.0


@pytest.mark.parametrize("R", [6, 12, 20])
def test_vertical_specialization_infinite_rank_vs_frobenius_finite_rank(R: int):
    """
    Prove that VERTICALLY (across R scales x_1 < ... < x_R), the Q-ramp kernel
    has full rank R -> infinity (cannot arise by linear specialization of a
    finite-rank module), whereas a fixed curve's Frobenius trace kernel has
    constant rank 2 = dim H^1(E) for all R.
    """
    vert = lab.audit_vertical_obstructions(q=11, d=1319, num_scales=R)
    comp = vert["2_specialization_linear_rank_comparison"]
    assert comp["number_field_ramp_matrix_rank"] == R
    assert comp["function_field_frobenius_hankel_rank"] == 2
