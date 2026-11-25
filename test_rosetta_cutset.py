"""
Pytest verification for `upstream_cutset_rosetta.py`, the live snapshot
`rosetta_cutset_snapshot.json` (`ce1560fb`), and Tab 6 of `docs/index.html`.
"""

from __future__ import annotations

import json
from pathlib import Path

import upstream_cutset_rosetta as ucr


SNAPSHOT_PATH = Path(__file__).resolve().parent / "rosetta_cutset_snapshot.json"


def _load_snapshot() -> dict:
    return json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))


def test_rosetta_stone_16_verified_permalinks_and_declarations():
    snap = _load_snapshot()
    assert snap["upstream_repo"] == "leanprover-community/riemann-hypothesis"
    assert snap["commit_short"] == "ce1560fb"
    assert len(snap["commit_sha"]) == 40

    rows = snap["rosetta_stone_16"]
    assert len(rows) == 16

    expected_labels = {
        "Main Headline",
        "Lemma 1.1 & 1.2",
        "Lemma 2.1",
        "(H0)",
        "(H1)",
        "(H2)",
        "(H3)",
        "(H4)",
        "(H5)",
        "(H6)",
        "(H7)",
        "(H8)",
        "(H9)",
        "(H10) & (H11)",
        "Estimates B & C",
        "(H12)",
    }
    assert {r["hypothesis"] for r in rows} == expected_labels

    by_hyp = {r["hypothesis"]: r for r in rows}
    assert by_hyp["Main Headline"]["line"] == 12
    assert "riemannHypothesis" in by_hyp["Main Headline"]["declaration_excerpt"]
    assert by_hyp["(H3)"]["line"] == 1214
    assert "structure LocHW" in by_hyp["(H3)"]["declaration_excerpt"]
    assert by_hyp["(H9)"]["line"] == 1230
    assert "vwin_of_row" in by_hyp["(H9)"]["declaration_excerpt"]

    for r in rows:
        assert r["line"] >= 1
        assert r["lean_path"].startswith("src/") and r["lean_path"].endswith(".lean")
        assert len(r["declaration_excerpt"].strip()) > 10
        expected_suffix = f"/{snap['commit_sha']}/{r['lean_path']}#L{r['line']}"
        assert r["github_permalink"].endswith(expected_suffix)


def test_abelian_variety_cutset_mycelial_mesh_invariants():
    """
    Verify the key architectural discovery inside `Development.AbelianVarietyTheory`:
    out of 25,765 reachable modules (6.14M lines), 10,292 unique modules (39.95%)
    are directly imported across 33,148 incoming cross-boundary edges from 5,653 callers.
    """
    snap = _load_snapshot()
    avt = snap["abelian_variety_cutset"]

    assert avt["total_avt_reachable_modules"] == 25765
    assert avt["total_avt_lines"] == 6143414
    assert avt["internal_avt_to_avt_edges"] == 134679
    assert avt["cross_boundary_incoming_edges"] == 33148
    assert avt["unique_avt_boundary_gateway_modules"] == 10292
    assert avt["gateway_compression_ratio_pct"] == 39.95
    assert avt["unique_external_caller_modules"] == 5653

    gateways = avt["top_20_boundary_gateways"]
    assert len(gateways) == 20
    assert gateways[0]["avt_module"] == "Development.AbelianVarietyTheory.Definitions.Divisors.MatchingOfParts"
    assert gateways[0]["external_callers_count"] == 247
    assert gateways[1]["avt_module"] == "Development.AbelianVarietyTheory.BlockTheory.DivisorBase.WeilTensorModels"
    assert gateways[1]["external_callers_count"] == 225


def test_inter_role_dependency_flows():
    snap = _load_snapshot()
    flows = snap["top_inter_role_dependency_flows"]
    assert len(flows) >= 10
    for f in flows:
        assert f["from_role"] != f["to_role"]
        assert f["import_edges"] > 0

    top = flows[0]
    assert top["from_role"] == "H12_Assembly_And_Windows"
    assert top["to_role"] == "H9_Galois_AbelianVarieties"
    assert top["import_edges"] == 23550


def test_rosetta_cutset_cli_snapshot_mode(capsys):
    rc = ucr.main(["--snapshot", str(SNAPSHOT_PATH)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "10,292 gateway modules (39.95% of AVT!)" in out
    assert "VwinOfRow.lean:L1230" in out


def test_github_pages_dashboard_contains_rosetta_panel():
    html = (Path(__file__).resolve().parent / "docs" / "index.html").read_text(encoding="utf-8")
    for required in (
        "panel-rosetta",
        "10,292",
        "39.95%",
        "VwinOfRow.lean",
        "MatchingOfParts",
        "ce1560fba43e59629abdf3c27a279218fad9f77a",
    ):
        assert required in html
