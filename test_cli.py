"""CLI integration tests for rh_audit.py."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

CLI = str(Path(__file__).resolve().parent / "rh_audit.py")


def _run_json(*args: str) -> dict:
    proc = subprocess.run(
        [sys.executable, CLI, *args, "--json"],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(proc.stdout)


def test_cli_spectral_double_zero():
    data = _run_json("spectral", "--num-zeros", "30", "--inject-double", "0")
    assert data["passes_one_sided_window_bound"] is True
    assert data["multiplicity_verdict"] == "MULTIPLE_ZERO_DETECTED"
    assert data["one_sided_blindness_triggered"] is True


def test_cli_curve_non_cm_and_cm():
    non_cm = _run_json("curve", "--limit", "2000")
    assert non_cm["hasse_weil_holds"] is True
    assert "NON_CM_GL2_SATO_TATE" in non_cm["distribution_verdict"]
    assert non_cm["gl1_table_captures_curve"] is False

    # CM curve y^2 = x^3 - x (j = 1728): half of primes are supersingular (a_p = 0)
    cm = _run_json("curve", "--a1", "0", "--a2", "0", "--a3", "0", "--a4", "-1", "--a6", "0", "--limit", "2000")
    assert cm["hasse_weil_holds"] is True
    assert "CM_ELLIPTIC_CURVE" in cm["distribution_verdict"]
    assert cm["supersingular_fraction"] > 0.45


def test_cli_table_d1319():
    data = _run_json("table", "--q", "11", "--d", "1319", "--N", "16")
    assert data["class_number_forms"] == 45
    assert data["identity_abs_error"] < 1e-10
    assert data["inequality_slack_over_log_h"] > 50.0


def test_cli_scale():
    data = _run_json("scale", "--x", "1e10", "--eps", "0.01")
    assert data["dominant_error_exponent"] == 0.500625
    assert data["certificate_dimensions_log10"]["log10_k_lattice_rank"] == 46.0
    assert data["fixed_object_comparison"]["elliptic_curve_dim_H1"] == 2


def test_cli_report():
    data = _run_json("report")
    assert set(data.keys()) == {
        "use_case_1_spectral_double_zero",
        "use_case_2_curve_E11a1",
        "use_case_3_table_d1319",
        "use_case_4_scale_1e10",
    }
