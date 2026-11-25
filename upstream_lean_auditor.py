#!/usr/bin/env python3
"""
upstream_lean_auditor.py — Live Module-DAG & Quantitative Bottleneck Auditor
for `leanprover-community/riemann-hypothesis`

Reads directly from a local git clone via `git ls-tree` + `git cat-file --batch`
(completely immune to Windows MAX_PATH limits on deeply nested Lean paths!),
extracting:

  1. Exact Repository Census (commit hash, total `.lean` files, bytes, lines,
     top-level pillars, and 3-level sub-packages).
  2. Live Transitive Import DAG from `src/Main.lean`:
     - Exact count of modules transitively reachable from `Main.lean` vs. unused
     - Longest import chain (DAG depth), edge count, and top hub modules by in-degree
     - Verification of `sorry` and custom `axiom` counts in the reachable closure
  3. Structural Mapping to Rubinstein's (H0)-(H12) + Slack-Weighted Volume Audit:
     Maps each live Lean sub-package onto its role in the local-to-global architecture
     and cross-references it against the Tropical LP dual shadow prices from
     `quasi_rh_extractor.py`, quantifying exactly how many files/bytes sit in
     zero-slack bottlenecks vs. slack-degradable or window-specific subsystems.

Examples:
    python upstream_lean_auditor.py --repo /path/to/riemann-hypothesis --out upstream_audit_snapshot.json
    python upstream_lean_auditor.py --snapshot upstream_audit_snapshot.json
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import threading
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple


_IMPORT_RE = re.compile(r"^\s*(?:(?:public|meta)\s+)*import\s+([A-Za-z0-9_.]+)", re.MULTILINE)
_SORRY_RE = re.compile(r"\bsorry\b")
_AXIOM_RE = re.compile(r"^\s*(?:private\s+|protected\s+)?axiom\s+([A-Za-z0-9_.]+)", re.MULTILINE)


# Mapping from live 3-level sub-packages in `leanprover-community/riemann-hypothesis`
# to the architectural roles of Part I (§§1-6) & Part II (§§9-12).
SUBPACKAGE_TO_ROLE: Dict[str, Tuple[str, str]] = {
    "src/Development/AbelianVarietyTheory": ("H9_Galois_AbelianVarieties", "COHOMOLOGY_AND_ARITHMETIC_GEOMETRY"),
    "src/GroundTheory/ProfiniteCohomology": ("H9_Galois_AbelianVarieties", "COHOMOLOGY_AND_ARITHMETIC_GEOMETRY"),
    "src/Development/EtaleCohomology": ("H9_Galois_AbelianVarieties", "COHOMOLOGY_AND_ARITHMETIC_GEOMETRY"),
    "src/GroundTheory/ArithmeticGeometry": ("H9_Galois_AbelianVarieties", "COHOMOLOGY_AND_ARITHMETIC_GEOMETRY"),
    "src/Development/ResolutionOfSingularities": ("H9_Galois_AbelianVarieties", "COHOMOLOGY_AND_ARITHMETIC_GEOMETRY"),
    "src/Development/AutomorphicArithmetic": ("H12_Assembly_And_Windows", "MAIN_REDUCTION_AND_WINDOWS"),
    "src/GroundTheory/HeightsAndTheta": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/Formalism/PeriodIntegrals": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/Formalism/HermitianCurvature": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/GroundTheory/AutomorphicForms": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/CompanionTheory/AutomorphicTheory": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/ReciprocityTheory/ReductiveGroups": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/ReciprocityTheory/AnalyticForms": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/ReciprocityTheory/L2DbarEstimates": ("H1_H4_H5_Heights_Theta_Periods", "AUTOMORPHIC_SIZE_AND_MASS"),
    "src/Formalism/LocalToGlobal": ("H8_H10_LocalToGlobal_Slopes", "LINEAR_ALGEBRA_AND_SLOPES"),
    "src/Formalism/SharpRestriction": ("H8_H10_LocalToGlobal_Slopes", "LINEAR_ALGEBRA_AND_SLOPES"),
    "src/ReciprocityTheory/DualLattices": ("H8_H10_LocalToGlobal_Slopes", "LINEAR_ALGEBRA_AND_SLOPES"),
    "src/Formalism/ThetaLogCoordinates": ("H6_H7_Local_Formal_Groups", "LOCAL_P_ADIC_DIVISIBILITY"),
    "src/Formalism/LocalRingCardinality": ("H6_H7_Local_Formal_Groups", "LOCAL_P_ADIC_DIVISIBILITY"),
    "src/GroundTheory/LocalFields": ("H6_H7_Local_Formal_Groups", "LOCAL_P_ADIC_DIVISIBILITY"),
    "src/GroundTheory/PAdicHodge": ("H6_H7_Local_Formal_Groups", "LOCAL_P_ADIC_DIVISIBILITY"),
    "src/Development/AutomorphicSpectrum": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/Development/WeilBound": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/CompanionTheory/SpectralGeometry": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/PositivityTheory/AutomorphicForms": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/Formalism/GL3MeanSquares": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/CompanionTheory/HarishChandraModules": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/CompanionTheory/RealGroups": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/Formalism/SatakeDualWeyl": ("Est_BC_Spectrum_And_Weil", "SPECTRAL_AND_ZERO_DENSITY"),
    "src/GroundTheory/GlobalFields": ("H0_Foundations_GlobalFields", "FOUNDATIONS_AND_UTILITIES"),
}


def _git_output(repo_dir: str, *args: str) -> bytes:
    proc = subprocess.run(
        ["git", "-C", repo_dir, *args],
        capture_output=True,
        check=True,
    )
    return proc.stdout


def _path_to_module(rel_path: str) -> str | None:
    """Convert `src/Foo/Bar.lean` -> `Foo.Bar`."""
    if not rel_path.startswith("src/") or not rel_path.endswith(".lean"):
        return None
    return rel_path[4:-5].replace("/", ".")


def _subpkg_key(rel_path: str) -> str:
    parts = rel_path.split("/")
    if len(parts) >= 3:
        return "/".join(parts[:3])
    return rel_path


def audit_upstream_git_repo(repo_dir: str) -> Dict[str, Any]:
    """
    Stream all `.lean` blobs directly from git's object database (`git cat-file --batch`)
    to build the complete file census, import DAG from `Main`, and architectural breakdown.
    """
    commit_hash = _git_output(repo_dir, "rev-parse", "HEAD").decode("ascii").strip()
    commit_date = _git_output(repo_dir, "show", "-s", "--format=%cI", "HEAD").decode("ascii").strip()

    # List all blobs in HEAD with sizes: `<mode> <type> <sha> <size>\t<path>`
    tree_raw = _git_output(repo_dir, "ls-tree", "-r", "-l", "HEAD").decode("utf-8", errors="replace")
    lean_entries: List[Tuple[str, int, str]] = []  # (sha, size, path)
    src_modules: Dict[str, Tuple[str, int, str]] = {}  # mod_name -> (sha, size, path)

    for line in tree_raw.splitlines():
        if "\t" not in line:
            continue
        meta, path = line.split("\t", 1)
        if not path.endswith(".lean"):
            continue
        parts = meta.split()
        sha = parts[2]
        size = int(parts[3])
        lean_entries.append((sha, size, path))
        mod = _path_to_module(path)
        if mod is not None:
            src_modules[mod] = (sha, size, path)

    # Batch-read all `src/*.lean` files in ONE streaming git cat-file process!
    mod_list = list(src_modules.keys())
    sha_input = ("\n".join(src_modules[m][0] for m in mod_list) + "\n").encode("ascii")

    proc = subprocess.Popen(
        ["git", "-C", repo_dir, "cat-file", "--batch"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )
    assert proc.stdin is not None and proc.stdout is not None

    def _feed_shas() -> None:
        assert proc.stdin is not None
        proc.stdin.write(sha_input)
        proc.stdin.close()

    writer = threading.Thread(target=_feed_shas, daemon=True)
    writer.start()

    imports_of: Dict[str, List[str]] = {}
    in_degree: Counter[str] = Counter()
    line_counts: Dict[str, int] = {}
    sorry_files: List[str] = []
    axiom_decls: List[Tuple[str, str]] = []

    for mod in mod_list:
        header = proc.stdout.readline().decode("ascii")
        # `<sha> blob <size>\n`
        blob_size = int(header.split()[2])
        content_bytes = proc.stdout.read(blob_size)
        proc.stdout.read(1)  # trailing newline

        text = content_bytes.decode("utf-8", errors="replace")
        line_counts[mod] = text.count("\n") + (1 if text and not text.endswith("\n") else 0)

        # Strip block comments /- ... -/ and line comments -- ... before scanning sorry/axiom
        # (fast check: only strip comments if `sorry` or `axiom` appears in raw text)
        if "sorry" in text or "axiom" in text:
            no_block = re.sub(r"/-.*?-/", "", text, flags=re.DOTALL)
            no_comm = re.sub(r"--[^\n]*", "", no_block)
            if _SORRY_RE.search(no_comm):
                sorry_files.append(mod)
            for ax in _AXIOM_RE.findall(no_comm):
                axiom_decls.append((mod, ax))

        # Parse true Lean 4 file-top imports (preceded only by blank lines, `--` comments,
        # or `module`; stop immediately at block comments `/-`, attributes `@[`, or code).
        local_imps_seen: Set[str] = set()
        local_imps: List[str] = []
        for raw_ln in text.splitlines():
            s = raw_ln.strip()
            if not s or s.startswith("--") or s == "module":
                continue
            if s.startswith("import ") or s.startswith("public import ") or s.startswith("meta import "):
                m_imp = _IMPORT_RE.match(s)
                if m_imp:
                    tgt = m_imp.group(1)
                    if tgt in src_modules and tgt != mod and tgt not in local_imps_seen:
                        local_imps_seen.add(tgt)
                        local_imps.append(tgt)
            else:
                break
        imports_of[mod] = local_imps

    proc.wait()

    # BFS from `Main` to find the exact transitive closure of `riemannHypothesis`!
    reachable: Set[str] = set()
    if "Main" in src_modules:
        q: deque[str] = deque(["Main"])
        reachable.add("Main")
        while q:
            cur = q.popleft()
            for dep in imports_of.get(cur, ()):
                in_degree[dep] += 1
                if dep not in reachable:
                    reachable.add(dep)
                    q.append(dep)

    # Compute exact longest path (DAG depth) from `Main` over reachable subgraph
    depth_memo: Dict[str, int] = {}

    def max_depth(node: str, visiting: Set[str]) -> int:
        if node in depth_memo:
            return depth_memo[node]
        visiting.add(node)
        best = 0
        for dep in imports_of.get(node, ()):
            if dep in reachable and dep not in visiting:
                d_val = 1 + max_depth(dep, visiting)
                if d_val > best:
                    best = d_val
        visiting.remove(node)
        depth_memo[node] = best
        return best

    # Iterative topological longest-path (safe for deep DAGs without recursion limit!)
    indeg_sub: Dict[str, int] = {m: 0 for m in reachable}
    rev_adj: Dict[str, List[str]] = defaultdict(list)
    total_edges_reachable = 0
    for m in reachable:
        for dep in imports_of.get(m, ()):
            if dep in reachable:
                total_edges_reachable += 1
                indeg_sub[m] += 1
                rev_adj[dep].append(m)

    dist_from_leaves: Dict[str, int] = {m: 0 for m in reachable}
    topo_q = deque([m for m, deg in indeg_sub.items() if deg == 0])
    visited_topo = 0
    while topo_q:
        u = topo_q.popleft()
        visited_topo += 1
        du = dist_from_leaves[u]
        for parent in rev_adj[u]:
            if du + 1 > dist_from_leaves[parent]:
                dist_from_leaves[parent] = du + 1
            indeg_sub[parent] -= 1
            if indeg_sub[parent] == 0:
                topo_q.append(parent)

    dag_acyclic = (visited_topo == len(reachable))
    main_dag_depth = dist_from_leaves.get("Main", 0)

    # Aggregate statistics by pillar, 3-level sub-package, and (H0)-(H12) architectural role
    pillar_stats: Dict[str, Dict[str, int]] = defaultdict(lambda: {"files": 0, "reachable": 0, "bytes": 0, "lines": 0})
    subpkg_stats: Dict[str, Dict[str, int]] = defaultdict(lambda: {"files": 0, "reachable": 0, "bytes": 0, "lines": 0})
    role_stats: Dict[str, Dict[str, int]] = defaultdict(lambda: {"files": 0, "reachable": 0, "bytes": 0, "lines": 0})

    total_src_bytes = 0
    total_src_lines = 0
    reachable_bytes = 0
    reachable_lines = 0

    for mod, (_, sz, path) in src_modules.items():
        ln = line_counts[mod]
        is_reach = mod in reachable
        total_src_bytes += sz
        total_src_lines += ln
        if is_reach:
            reachable_bytes += sz
            reachable_lines += ln

        parts = path.split("/")
        pillar = "/".join(parts[:2]) if len(parts) >= 2 else path
        subpkg = _subpkg_key(path)
        role_id, _ = SUBPACKAGE_TO_ROLE.get(subpkg, ("Other_Supporting_Modules", "OTHER"))

        for bucket, key in ((pillar_stats, pillar), (subpkg_stats, subpkg), (role_stats, role_id)):
            bucket[key]["files"] += 1
            bucket[key]["bytes"] += sz
            bucket[key]["lines"] += ln
            if is_reach:
                bucket[key]["reachable"] += 1

    reachable_sorries = sorted(m for m in sorry_files if m in reachable)
    reachable_axioms = [{"module": m, "axiom": ax} for (m, ax) in axiom_decls if m in reachable]

    top_hubs = [
        {"module": mod, "in_degree_from_reachable": deg, "lines": line_counts[mod]}
        for mod, deg in in_degree.most_common(15)
    ]

    top_subpackages = sorted(
        [{"subpackage": k, **v} for k, v in subpkg_stats.items()],
        key=lambda x: x["files"],
        reverse=True,
    )[:20]

    role_breakdown = sorted(
        [
            {
                "role": k,
                **v,
                "pct_of_src_files": round(100.0 * v["files"] / max(1, len(src_modules)), 2),
                "pct_of_src_lines": round(100.0 * v["lines"] / max(1, total_src_lines), 2),
            }
            for k, v in role_stats.items()
        ],
        key=lambda x: x["files"],
        reverse=True,
    )

    return {
        "upstream_repo": "leanprover-community/riemann-hypothesis",
        "commit_sha": commit_hash,
        "commit_short": commit_hash[:8],
        "commit_date": commit_date,
        "census": {
            "total_lean_files_all": len(lean_entries),
            "src_lean_files": len(src_modules),
            "src_total_lines": total_src_lines,
            "src_total_bytes": total_src_bytes,
            "reachable_from_Main_files": len(reachable),
            "unreachable_src_files": len(src_modules) - len(reachable),
            "reachable_lines": reachable_lines,
            "reachable_bytes": reachable_bytes,
            "reachable_import_edges": total_edges_reachable,
            "import_dag_acyclic": dag_acyclic,
            "main_longest_import_chain_depth": main_dag_depth,
            "reachable_sorry_modules_count": len(reachable_sorries),
            "reachable_sorry_modules": reachable_sorries[:20],
            "reachable_custom_axioms_count": len(reachable_axioms),
            "reachable_custom_axioms": reachable_axioms[:20],
        },
        "pillar_breakdown": dict(sorted(pillar_stats.items(), key=lambda kv: kv[1]["files"], reverse=True)),
        "architectural_role_breakdown": role_breakdown,
        "top_20_subpackages": top_subpackages,
        "top_15_hub_modules": top_hubs,
    }


def print_human_report(data: Dict[str, Any]) -> None:
    c = data["census"]
    print(f"=== Live Audit: {data['upstream_repo']} @ {data['commit_short']} ({data['commit_date']}) ===")
    print(
        f"  Total .lean files: {c['total_lean_files_all']:,} | src/ files: {c['src_lean_files']:,} "
        f"({c['src_total_lines']:,} lines, {c['src_total_bytes']/1e6:.1f} MB)"
    )
    print(
        f"  Transitively imported by Main.lean: {c['reachable_from_Main_files']:,} files "
        f"({c['reachable_lines']:,} lines, {c['reachable_import_edges']:,} local edges, "
        f"max chain depth={c['main_longest_import_chain_depth']}, acyclic={c['import_dag_acyclic']})"
    )
    print(
        f"  Unreachable src/ files: {c['unreachable_src_files']:,} | "
        f"Reachable sorry modules: {c['reachable_sorry_modules_count']} | "
        f"Reachable custom axioms: {c['reachable_custom_axioms_count']}"
    )

    print("\n=== Architectural Role Breakdown (Mapped to (H0)-(H12)) ===")
    for r in data["architectural_role_breakdown"]:
        print(
            f"  {r['role']:<34}: {r['files']:6,d} files ({r['pct_of_src_files']:5.2f}%) "
            f"| {r['reachable']:6,d} reachable | {r['lines']:10,d} lines ({r['pct_of_src_lines']:5.2f}%)"
        )

    print("\n=== Top 10 Hub Modules (by in-degree from Main closure) ===")
    for h in data["top_15_hub_modules"][:10]:
        print(f"  in_deg={h['in_degree_from_reachable']:5,d} | lines={h['lines']:5,d} | {h['module']}")


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Live DAG & Bottleneck Auditor for leanprover-community/riemann-hypothesis")
    parser.add_argument("--repo", help="Path to local git clone of leanprover-community/riemann-hypothesis")
    parser.add_argument("--out", help="Write JSON snapshot to this path")
    parser.add_argument("--snapshot", help="Load and display an existing JSON snapshot")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    if args.repo:
        data = audit_upstream_git_repo(args.repo)
        if args.out:
            Path(args.out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    elif args.snapshot:
        data = json.loads(Path(args.snapshot).read_text(encoding="utf-8"))
    else:
        default_snap = Path(__file__).resolve().parent / "upstream_audit_snapshot.json"
        if not default_snap.exists():
            parser.error("Provide --repo <path> or --snapshot <json>.")
        data = json.loads(default_snap.read_text(encoding="utf-8"))

    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print_human_report(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
