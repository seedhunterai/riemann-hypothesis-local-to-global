"""
Pytest verification for `upstream_lean_auditor.py` and the live snapshot
`upstream_audit_snapshot.json` of `leanprover-community/riemann-hypothesis` @ ce1560fb.
"""

from __future__ import annotations

import json
from pathlib import Path

import upstream_lean_auditor as ua


SNAPSHOT_PATH = Path(__file__).resolve().parent / "upstream_audit_snapshot.json"


def _load_snapshot() -> dict:
    return json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))


def test_upstream_snapshot_integrity_and_acyclic_dag():
    snap = _load_snapshot()
    c = snap["census"]

    assert snap["upstream_repo"] == "leanprover-community/riemann-hypothesis"
    assert snap["commit_short"] == "ce1560fb"
    assert c["total_lean_files_all"] == 58880
    assert c["src_lean_files"] == 58871
    assert c["reachable_from_Main_files"] == 58858
    assert c["reachable_from_Main_files"] + c["unreachable_src_files"] == c["src_lean_files"]
    assert c["import_dag_acyclic"] is True
    assert c["main_longest_import_chain_depth"] == 95
    assert c["reachable_sorry_modules_count"] == 0
    assert c["reachable_custom_axioms_count"] == 0


def test_upstream_architectural_role_accounting_and_asymmetry():
    """
    Verify exact partition accounting across all 58,871 src/ files, and confirm
    the massive architectural asymmetry:
      - H9 (Abelian varieties & Galois/etale/profinite cohomology) > 45% of the proof
      - H6-H7 (Local Lubin-Tate / formal groups) < 1% of the proof!
    """
    snap = _load_snapshot()
    c = snap["census"]
    roles = {r["role"]: r for r in snap["architectural_role_breakdown"]}

    assert sum(r["files"] for r in roles.values()) == c["src_lean_files"]
    assert sum(r["lines"] for r in roles.values()) == c["src_total_lines"]
    assert sum(r["reachable"] for r in roles.values()) == c["reachable_from_Main_files"]

    assert roles["H9_Galois_AbelianVarieties"]["pct_of_src_files"] > 45.0
    assert roles["H6_H7_Local_Formal_Groups"]["pct_of_src_files"] < 1.0


def test_upstream_top_hub_modules_present():
    snap = _load_snapshot()
    hub_names = {h["module"] for h in snap["top_15_hub_modules"]}
    assert "Formalism.LocalToGlobal.IndefiniteReal.LinksDefs" in hub_names
    assert "Formalism.LocalToGlobal.ProjectiveLattices.ErealAddIInf" in hub_names
    assert "Formalism.SharpRestriction.GlobalSchurBounds.GeneralLinearActions.GenusMassVolumes" in hub_names


def test_upstream_cli_snapshot_mode(capsys):
    rc = ua.main(["--snapshot", str(SNAPSHOT_PATH)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "ce1560fb" in out
    assert "max chain depth=95, acyclic=True" in out


def test_github_pages_dashboard_contains_all_five_panels():
    html = (Path(__file__).resolve().parent / "docs" / "index.html").read_text(encoding="utf-8")
    for required in ("panel-upstream", "panel-quasirh", "panel-fitting", "panel-tower", "panel-zd", "ce1560fb", "58,880"):
        assert required in html
