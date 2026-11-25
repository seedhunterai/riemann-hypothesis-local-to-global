"""
Pytest verification for the frontier GL_2 local-to-global solver:
  - Exact analytic rank & leading central L-derivative L^(r)(E, 1) from local a_p
  - Tight Weil-Mestre explicit-formula GRH rank ceiling for ranks 0, 1, 2, and 3
"""

from __future__ import annotations

import pytest

import gl2_local_to_global as gl2


# Reference leading central derivatives L^(r)(E, 1) from LMFDB / Cremona tables:
REFERENCE_CURVES = [
    ("11a1", 0, 0.25384186),
    ("37a1", 1, 0.30599977),
    ("389a1", 2, 1.51863300),
    ("5077a1", 3, 10.39109940),
]


@pytest.mark.parametrize("label,expected_rank,expected_L_val", REFERENCE_CURVES)
def test_gl2_local_to_global_exact_rank_and_weil_bound(
    label: str, expected_rank: int, expected_L_val: float
):
    curve = gl2.BENCHMARK_CURVES[label]
    res = gl2.analyze_curve(curve, delta=7.5)

    exact = res["exact_local_to_global_solver"]
    weil = res["weil_mestre_explicit_formula"]

    # 1. Exact analytic rank solved purely from local Frobenius traces a_p
    assert exact["analytic_rank"] == expected_rank
    assert abs(exact["leading_L_derivative"] - expected_L_val) < 1e-4

    # 2. Weil-Mestre explicit formula tightly certifies r <= expected_rank (no slack!)
    assert weil["certified_integer_rank_ceiling"] == expected_rank
    assert res["rank_pinpointed_exactly"] is True
