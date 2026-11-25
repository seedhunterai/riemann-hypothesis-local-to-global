#!/usr/bin/env python3
"""
upstream_dominator_and_challenge.py — Frontier #1 & #2:
  1. Exact Dominator-Tree PR Target Finder & AGENTS.md Code-Health Census
     on `leanprover-community/riemann-hypothesis` (`ce1560fb`, 58,858 reachable modules).
     Finds mathematically certified 1-Vertex-Cut Encapsulated Islands (`side_entrances = 0`):
     subtrees `D(d)` where every path from `Main` into `D(d)` passes through `d`,
     so replacing/refactoring `d` unlocks deleting all `|D(d)| - 1` internal files
     in one atomic PR with zero broken imports elsewhere.
  2. `src/Challenge.lean` / `src/Main.lean` 5-Theorem Dependency Cone Split
     (`riemannHypothesis` vs. `DirichletLNonvanishingRightOfHalf` & the `HalfPlane.lean`
     Siegel-Weil detour) + Exact Q-Farkas Conductor-Aspect (`q = x^{alpha_q}`) Quasi-GRH Solver.

Examples:
    python upstream_dominator_and_challenge.py --repo /path/to/riemann-hypothesis --out dominator_challenge_snapshot.json
    python upstream_dominator_and_challenge.py --snapshot dominator_challenge_snapshot.json
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import threading
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import numpy as np
from scipy.optimize import linprog

import upstream_lean_auditor as ua


NOV20_BASELINE_FILES = 175_187
NOV20_BASELINE_LINES = 17_374_215

_DECL_HEAD_RE = re.compile(
    r"^\s*(?:public\s+|noncomputable\s+|protected\s+)*(theorem|lemma|def|structure)\s+([A-Za-z0-9_'.]+)"
)


def solve_conductor_quasi_grh(tau_k: float | Fraction, alpha_q: float | Fraction) -> Dict[str, Any]:
    """
    Exact Q-Farkas + Primal LP solver for the Conductor-Augmented Window-Free
    Quasi-GRH Blueprint Model (`q = x^{alpha_q}`, `alpha_q >= 0`, `tau_k >= 3/2`):

      Constraints on `(theta, alpha_m, alpha_k, alpha_n) >= 0`:
        1. end_chain       : 1 - alpha_m <= theta
        2. step_accum      : alpha_m <= theta
        3. rank_spread     : 1 + alpha_m - alpha_k <= theta
        4. table_collision : alpha_k - 2*alpha_n <= theta
        5. codim_span_q    : 1 + alpha_q + 2*alpha_n - alpha_k <= theta
        6. table_modulus_q : 5/2 + alpha_q - alpha_n <= theta
        7. twist_cap       : alpha_k <= tau_k

      Exact Closed-Form Theorem (proved in Lean 4 `Theorems.lean`):
        theta*(tau_k, alpha_q) = max((1 + alpha_q)/2,  2 + alpha_q - tau_k/3)
    """
    tk = Fraction(str(tau_k)) if not isinstance(tau_k, Fraction) else tau_k
    aq = Fraction(str(alpha_q)) if not isinstance(alpha_q, Fraction) else alpha_q
    if aq < 0:
        raise ValueError("alpha_q must be >= 0")
    if tk < Fraction(3, 2):
        raise ValueError("tau_k must be >= 3/2 for the two-branch closed form")

    branch_conv = Fraction(1, 2) * (1 + aq)
    branch_twist = Fraction(2, 1) + aq - tk / 3
    theta_exact = max(branch_conv, branch_twist)

    # Exact primal witness and exact Farkas dual multipliers
    kink = Fraction(9, 2) + Fraction(3, 2) * aq
    if tk <= kink:
        am = Fraction(1, 2)
        ak = tk
        an = Fraction(1, 2) + tk / 3
        duals = {
            "end_chain": "0",
            "step_accum": "0",
            "rank_spread": "0",
            "table_collision": "0",
            "codim_span_q": "1/3",
            "table_modulus_q": "2/3",
            "twist_cap": "1/3",
        }
    else:
        am = Fraction(1, 2)
        ak = kink
        an = Fraction(2, 1) + aq / 2
        duals = {
            "end_chain": "0",
            "step_accum": "0",
            "rank_spread": "0",
            "table_collision": "1/2",
            "codim_span_q": "1/2",
            "table_modulus_q": "0",
            "twist_cap": "0",
        }

    # Verify all 7 inequalities in exact Fraction arithmetic
    assert 1 - am <= theta_exact
    assert am <= theta_exact
    assert 1 + am - ak <= theta_exact
    assert ak - 2 * an <= theta_exact
    assert 1 + aq + 2 * an - ak <= theta_exact
    assert Fraction(5, 2) + aq - an <= theta_exact
    assert ak <= tk

    # Cross-check against numerical HiGHS LP
    tk_f, aq_f = float(tk), float(aq)
    # Variables: [theta, am, ak, an] >= 0
    A_ub = np.array([
        [-1.0, -1.0,  0.0,  0.0],  # -theta - am <= -1
        [-1.0,  1.0,  0.0,  0.0],  # -theta + am <= 0
        [-1.0,  1.0, -1.0,  0.0],  # -theta + am - ak <= -1
        [-1.0,  0.0,  1.0, -2.0],  # -theta + ak - 2*an <= 0
        [-1.0,  0.0, -1.0,  2.0],  # -theta - ak + 2*an <= -1 - aq
        [-1.0,  0.0,  0.0, -1.0],  # -theta - an <= -2.5 - aq
        [ 0.0,  0.0,  1.0,  0.0],  # ak <= tk
    ])
    b_ub = np.array([-1.0, 0.0, -1.0, 0.0, -1.0 - aq_f, -2.5 - aq_f, tk_f])
    lp = linprog(c=[1.0, 0.0, 0.0, 0.0], A_ub=A_ub, b_ub=b_ub, bounds=[(0.0, None)] * 4, method="highs")
    lp_theta = float(lp.x[0])
    assert abs(lp_theta - float(theta_exact)) < 1e-9

    if aq >= 1:
        regime = "CONDUCTOR_CONVEXITY_BARRIER (alpha_q >= 1 => theta* >= 1)"
    elif theta_exact >= 1:
        regime = "INSUFFICIENT_TWIST_CAP (theta* >= 1)"
    elif theta_exact == Fraction(1, 2):
        regime = "FULL_GRH_AT_TRIVIAL_CONDUCTOR"
    else:
        regime = "NONTRIVIAL_WINDOW_FREE_QUASI_GRH (theta* < 1)"

    return {
        "tau_k": str(tk),
        "alpha_q": str(aq),
        "theta_star_exact": str(theta_exact),
        "theta_star_float": round(float(theta_exact), 6),
        "lp_theta_float": round(lp_theta, 6),
        "duality_gap": float(f"{abs(lp_theta - float(theta_exact)):.2e}"),
        "nontrivial_power_saving": bool(theta_exact < 1),
        "regime": regime,
        "primal_witness_exact": {
            "alpha_m": str(am),
            "alpha_k": str(ak),
            "alpha_n": str(an),
        },
        "farkas_duals_exact": duals,
    }


def extract_dominators_and_challenge(repo_dir: str) -> Dict[str, Any]:
    commit_hash = ua._git_output(repo_dir, "rev-parse", "HEAD").decode("ascii").strip()
    short_sha = commit_hash[:8]

    tree_raw = ua._git_output(repo_dir, "ls-tree", "-r", "-l", "HEAD").decode("utf-8", errors="replace")
    src_modules: Dict[str, Tuple[str, int, str]] = {}
    total_lean_all = 0
    for line in tree_raw.splitlines():
        if "\t" not in line:
            continue
        meta, path = line.split("\t", 1)
        if path.endswith(".lean"):
            total_lean_all += 1
        mod = ua._path_to_module(path)
        if mod is not None:
            parts = meta.split()
            src_modules[mod] = (parts[2], int(parts[3]), path)

    mod_list = list(src_modules.keys())
    sha_input = ("\n".join(src_modules[m][0] for m in mod_list) + "\n").encode("ascii")

    proc = subprocess.Popen(
        ["git", "-C", repo_dir, "cat-file", "--batch"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )
    assert proc.stdin is not None and proc.stdout is not None

    def _feed() -> None:
        assert proc.stdin is not None
        proc.stdin.write(sha_input)
        proc.stdin.close()

    threading.Thread(target=_feed, daemon=True).start()

    imports_of: Dict[str, List[str]] = {}
    line_counts: Dict[str, int] = {}
    decls_last: Dict[str, List[str]] = {}

    clump_member_markers = 0
    clumped_files_count = 0
    backward_opt_files = 0
    exact_by_kernel_files = 0

    for mod in mod_list:
        header = proc.stdout.readline().decode("ascii")
        blob_size = int(header.split()[2])
        content = proc.stdout.read(blob_size).decode("utf-8", errors="replace")
        proc.stdout.read(1)

        line_counts[mod] = content.count("\n") + (1 if content and not content.endswith("\n") else 0)

        # AGENTS.md code-health & clumping telemetry
        cm = content.count("end clump member:")
        if cm > 0:
            clumped_files_count += 1
            clump_member_markers += cm
        if "set_option backward" in content:
            backward_opt_files += 1
        if "exact_by_kernel" in content:
            exact_by_kernel_files += 1

        seen: Set[str] = set()
        imps: List[str] = []
        for raw_ln in content.splitlines():
            s = raw_ln.strip()
            if not s or s.startswith("--") or s == "module":
                continue
            if s.startswith(("import ", "public import ", "meta import ")):
                m_imp = ua._IMPORT_RE.match(s)
                if m_imp:
                    tgt = m_imp.group(1)
                    if tgt in src_modules and tgt != mod and tgt not in seen:
                        seen.add(tgt)
                        imps.append(tgt)
            else:
                break
        imports_of[mod] = imps

        d_names = []
        for raw_ln in content.splitlines():
            m_d = _DECL_HEAD_RE.match(raw_ln)
            if m_d:
                d_names.append(f"{m_d.group(1)} {m_d.group(2)}")
        decls_last[mod] = d_names[-3:]

    proc.wait()

    def cone_from(roots: List[str]) -> Set[str]:
        vis: Set[str] = set()
        q: deque[str] = deque()
        for r in roots:
            if r in src_modules and r not in vis:
                vis.add(r)
                q.append(r)
        while q:
            u = q.popleft()
            for v in imports_of.get(u, ()):
                if v not in vis:
                    vis.add(v)
                    q.append(v)
        return vis

    # --- Part A: Challenge.lean / Main.lean 5-Theorem Dependency Cone Split ---
    reach_main = cone_from(["Main"])
    cone_statements = cone_from(["Statements"])
    cone_claim_strip = cone_from(["Development.AutomorphicArithmetic.CharacterBoxes.ClaimStrip"])
    cone_half_plane = cone_from(["Development.AutomorphicArithmetic.CharacterBoxes.HalfPlane"])
    cone_meromorphic_detour = cone_from([
        "Development.WeilBound.CoreStatements.AutomorphicKernels.SiegelWeilFormula.ZetaEulerProducts.MeromorphicAtProd"
    ])
    cone_main_conclusion = cone_from(["Development.AutomorphicArithmetic.CharacterBoxes.MainConclusion"])

    strip_exclusive_outside_main_conc = cone_claim_strip - cone_main_conclusion
    detour_overlap_with_main_conc = cone_meromorphic_detour & cone_main_conclusion

    # --- Part B: Exact Dominator Tree on the Acyclic Reachable DAG from `Main` ---
    # Topological sort on induced subgraph `reach_main`
    indeg: Dict[str, int] = {m: 0 for m in reach_main}
    parents: Dict[str, List[str]] = defaultdict(list)
    for u in reach_main:
        for v in imports_of.get(u, ()):
            if v in reach_main:
                indeg[v] += 1
                parents[v].append(u)

    topo: List[str] = []
    dq: deque[str] = deque([m for m, d in indeg.items() if d == 0])
    while dq:
        u = dq.popleft()
        topo.append(u)
        for v in imports_of.get(u, ()):
            if v in reach_main:
                indeg[v] -= 1
                if indeg[v] == 0:
                    dq.append(v)

    assert len(topo) == len(reach_main), f"Expected DAG: topo={len(topo)} vs reach={len(reach_main)}"

    idom: Dict[str, str] = {"Main": "Main"}
    dom_depth: Dict[str, int] = {"Main": 0}

    def lca(a: str, b: str) -> str:
        while dom_depth[a] > dom_depth[b]:
            a = idom[a]
        while dom_depth[b] > dom_depth[a]:
            b = idom[b]
        while a != b:
            a = idom[a]
            b = idom[b]
        return a

    for v in topo:
        if v == "Main":
            continue
        ps = parents[v]
        cur_dom = ps[0]
        for p in ps[1:]:
            cur_dom = lca(cur_dom, p)
        idom[v] = cur_dom
        dom_depth[v] = dom_depth[cur_dom] + 1

    # Build children in Dominator Tree & accumulate exact subtree sizes bottom-up
    dom_children: Dict[str, List[str]] = defaultdict(list)
    for v, d in idom.items():
        if v != "Main":
            dom_children[d].append(v)

    dom_size: Dict[str, int] = {m: 1 for m in reach_main}
    dom_lines: Dict[str, int] = {m: line_counts[m] for m in reach_main}
    for v in reversed(topo):
        if v != "Main":
            p = idom[v]
            dom_size[p] += dom_size[v]
            dom_lines[p] += dom_lines[v]

    def get_role(m: str) -> str:
        _, _, p = src_modules[m]
        sp = ua._subpkg_key(p)
        return ua.SUBPACKAGE_TO_ROLE.get(sp, ("Other_Supporting_Modules", "OTHER"))[0]

    def collect_dom_subtree(root_d: str) -> Set[str]:
        s: Set[str] = {root_d}
        st = [root_d]
        while st:
            x = st.pop()
            for ch in dom_children.get(x, ()):
                s.add(ch)
                st.append(ch)
        return s

    def build_island_record(d: str) -> Dict[str, Any]:
        sub = collect_dom_subtree(d)
        # Verify 0 side entrances into `sub \ {d}` and count external foundation deps
        side_entrances = 0
        incoming_to_gateway = 0
        for p in parents[d]:
            if p not in sub:
                incoming_to_gateway += 1
        for w in sub:
            if w == d:
                continue
            for p in parents[w]:
                if p not in sub:
                    side_entrances += 1

        ext_deps: Set[str] = set()
        for w in sub:
            for dep in imports_of.get(w, ()):
                if dep in reach_main and dep not in sub:
                    ext_deps.add(dep)

        assert side_entrances == 0, f"Dominator invariant violated at {d}!"
        _, _, rel_p = src_modules[d]
        return {
            "gateway_module": d,
            "path": rel_p,
            "role": get_role(d),
            "dominated_files": dom_size[d],
            "deletable_internal_files": dom_size[d] - 1,
            "dominated_lines": dom_lines[d],
            "gateway_lines": line_counts[d],
            "incoming_callers_to_gateway": incoming_to_gateway,
            "side_entrances_to_interior": side_entrances,
            "external_foundation_deps": len(ext_deps),
            "headline_declarations": decls_last.get(d, []),
            "github_url": f"https://github.com/leanprover-community/riemann-hypothesis/blob/{commit_hash}/{rel_p}",
        }

    # Select maximal non-overlapping PR-sized islands (25 <= dominated_files <= 600):
    # A candidate `d` is "maximal PR-sized" if `25 <= dom_size[d] <= 600` and `dom_size[idom[d]] > 600`.
    pr_candidates = [
        m for m in reach_main
        if m != "Main" and 25 <= dom_size[m] <= 600 and dom_size[idom[m]] > 600
    ]
    pr_candidates.sort(key=lambda m: dom_size[m], reverse=True)
    top_pr_islands = [build_island_record(m) for m in pr_candidates[:20]]

    # Also find the top maximal encapsulated islands strictly inside `AbelianVarietyTheory` (`H9`)
    avt_candidates = [
        m for m in reach_main
        if m.startswith("Development.AbelianVarietyTheory.")
        and dom_size[m] >= 15
        and (not idom[m].startswith("Development.AbelianVarietyTheory.") or dom_size[idom[m]] > 500)
    ]
    avt_candidates.sort(key=lambda m: dom_size[m], reverse=True)
    top_avt_islands = [build_island_record(m) for m in avt_candidates[:12]]

    # Also record the Macro Dominator Spine (`dom_size >= 1000`)
    macro_dominators = sorted(
        [m for m in reach_main if dom_size[m] >= 1000],
        key=lambda m: dom_size[m],
        reverse=True,
    )
    macro_spine = [
        {
            "module": m,
            "idom": idom[m],
            "dominated_files": dom_size[m],
            "dominated_lines": dom_lines[m],
            "pct_of_reachable": round(100.0 * dom_size[m] / len(reach_main), 2),
        }
        for m in macro_dominators
    ]

    # Estimated pre-clump module count: unclumped files + total clump members
    unclumped_src_files = len(src_modules) - clumped_files_count
    reconstructed_preclump_modules = unclumped_src_files + clump_member_markers

    # Sample grid for Conductor-Aspect Quasi-GRH
    q_grh_grid = [
        solve_conductor_quasi_grh(tau_k=tk, alpha_q=aq)
        for aq in ("0", "1/4", "1/2", "3/4", "1")
        for tk in ("3", "18/5", "9/2", "6")
    ]

    return {
        "upstream_repo": "leanprover-community/riemann-hypothesis",
        "commit_sha": commit_hash,
        "commit_short": short_sha,
        "shrinkage_and_code_health": {
            "nov20_original_files": NOV20_BASELINE_FILES,
            "nov20_original_lines": NOV20_BASELINE_LINES,
            "current_all_lean_files": total_lean_all,
            "current_src_files": len(src_modules),
            "current_src_lines": sum(line_counts.values()),
            "file_reduction_pct_vs_nov20": round(
                100.0 * (NOV20_BASELINE_FILES - total_lean_all) / NOV20_BASELINE_FILES, 2
            ),
            "line_reduction_pct_vs_nov20": round(
                100.0 * (NOV20_BASELINE_LINES - sum(line_counts.values())) / NOV20_BASELINE_LINES, 2
            ),
            "clumped_files_count": clumped_files_count,
            "clump_member_markers_count": clump_member_markers,
            "reconstructed_preclump_src_modules": reconstructed_preclump_modules,
            "agents_md_debt_set_option_backward_files": backward_opt_files,
            "agents_md_debt_exact_by_kernel_files": exact_by_kernel_files,
        },
        "challenge_5_theorem_split": {
            "comparator_targets": [
                {
                    "theorem": "LFunction_ne_zero_of_half_lt_re",
                    "statement_type": "DirichletLNonvanishingRightOfHalf",
                    "proof_term_in_Main": "UpperTheory.FinalTheorems.everyCharacter",
                    "direct_provider": "Development.AutomorphicArithmetic.CharacterBoxes.MainConclusion",
                    "reachable_src_files": len(cone_main_conclusion),
                    "reachable_lines": sum(line_counts[m] for m in cone_main_conclusion),
                },
                {
                    "theorem": "riemannHypothesis",
                    "statement_type": "RiemannHypothesis",
                    "proof_term_in_Main": "Claim.zeta LFunction_ne_zero_of_half_lt_re",
                    "direct_provider": "ClaimStrip + MainConclusion (specializes N=1, chi=1!)",
                    "reachable_src_files": len(cone_main_conclusion | cone_claim_strip),
                    "reachable_lines": sum(line_counts[m] for m in (cone_main_conclusion | cone_claim_strip)),
                },
                {
                    "theorem": "riemannZeta_eq_zero_re",
                    "statement_type": "RiemannZetaZerosOnLine",
                    "proof_term_in_Main": "riemannHypothesis",
                    "direct_provider": "Definitional alias of riemannHypothesis",
                    "reachable_src_files": len(cone_main_conclusion | cone_claim_strip),
                    "reachable_lines": sum(line_counts[m] for m in (cone_main_conclusion | cone_claim_strip)),
                },
                {
                    "theorem": "LFunction_eq_zero_re",
                    "statement_type": "DirichletLZerosOnLine",
                    "proof_term_in_Main": "Claim.strip LFunction_ne_zero_of_half_lt_re",
                    "direct_provider": "ClaimStrip + MainConclusion (functional-equation reflection)",
                    "reachable_src_files": len(cone_main_conclusion | cone_claim_strip),
                    "reachable_lines": sum(line_counts[m] for m in (cone_main_conclusion | cone_claim_strip)),
                },
                {
                    "theorem": "LFunction_eq_zero_re_of_isPrimitive",
                    "statement_type": "PrimitiveDirichletLZerosOnLine",
                    "proof_term_in_Main": "exact LFunction_eq_zero_re N chi s h0 h1 hz",
                    "direct_provider": "2-line specialization of LFunction_eq_zero_re",
                    "reachable_src_files": len(cone_main_conclusion | cone_claim_strip),
                    "reachable_lines": sum(line_counts[m] for m in (cone_main_conclusion | cone_claim_strip)),
                },
            ],
            "statements_cone_files": len(cone_statements),
            "claim_strip_total_cone_files": len(cone_claim_strip),
            "claim_strip_total_cone_lines": sum(line_counts[m] for m in cone_claim_strip),
            "half_plane_meromorphic_detour_cone_files": len(cone_meromorphic_detour),
            "half_plane_meromorphic_detour_cone_lines": sum(line_counts[m] for m in cone_meromorphic_detour),
            "claim_strip_exclusive_outside_main_conclusion_files": sorted(strip_exclusive_outside_main_conc),
            "detour_overlap_with_main_conclusion_files": len(detour_overlap_with_main_conc),
        },
        "dominator_pr_finder": {
            "total_reachable_nodes": len(reach_main),
            "max_dominator_tree_depth": max(dom_depth.values()),
            "nontrivial_dominators_ge_2_files": sum(1 for v in dom_size.values() if v >= 2),
            "pr_sized_islands_25_to_600_count": len(pr_candidates),
            "macro_dominator_spine": macro_spine,
            "top_20_pr_ready_islands": top_pr_islands,
            "top_12_avt_h9_encapsulated_islands": top_avt_islands,
        },
        "conductor_aspect_quasi_grh": {
            "closed_form_formula": "theta*(tau_k, alpha_q) = max((1 + alpha_q)/2, 2 + alpha_q - tau_k/3)",
            "kink_threshold": "tau_k = (9 + 3*alpha_q)/2",
            "nontrivial_condition": "alpha_q < 1 and tau_k > 3 + 3*alpha_q",
            "grid_certificates": q_grh_grid,
        },
    }


def print_report(data: Dict[str, Any]) -> None:
    sh = data["shrinkage_and_code_health"]
    ch = data["challenge_5_theorem_split"]
    dm = data["dominator_pr_finder"]

    print(f"=== 1. Shrinkage & AGENTS.md Code-Health Census ({data['commit_short']}) ===")
    print(
        f"  Nov 20 Baseline -> Current : {sh['nov20_original_files']:,} -> {sh['current_all_lean_files']:,} files "
        f"(-{sh['file_reduction_pct_vs_nov20']}%) | {sh['nov20_original_lines']:,} -> {sh['current_src_lines']:,} lines (-{sh['line_reduction_pct_vs_nov20']}%)"
    )
    print(
        f"  Auto-Clumping Telemetry    : {sh['clumped_files_count']:,} clumped files holding {sh['clump_member_markers_count']:,} clump members "
        f"(reconstructed pre-clump pool: {sh['reconstructed_preclump_src_modules']:,} modules)"
    )
    print(
        f"  AGENTS.md Debt Markers     : set_option backward={sh['agents_md_debt_set_option_backward_files']:,} files | "
        f"exact_by_kernel={sh['agents_md_debt_exact_by_kernel_files']:,} files"
    )

    print("\n=== 2. Challenge.lean / Main.lean 5-Theorem Dependency Cone Split ===")
    print(
        f"  MainConclusion cone        : {ch['comparator_targets'][0]['reachable_src_files']:,} files "
        f"(proves DirichletLNonvanishingRightOfHalf FIRST; riemannHypothesis specializes N=1, chi=1!)"
    )
    print(
        f"  ClaimStrip cone            : {ch['claim_strip_total_cone_files']:,} files ({ch['claim_strip_total_cone_lines']:,} lines) — "
        f"HalfPlane.lean drags {ch['half_plane_meromorphic_detour_cone_files']:,} Siegel-Weil files via MeromorphicAtProd!"
    )

    print("\n=== 3. Exact Dominator-Tree PR Target Finder (0 Side-Entrances Certified) ===")
    print(
        f"  Dominator tree stats       : {dm['nontrivial_dominators_ge_2_files']:,} multi-file islands | "
        f"max dom depth={dm['max_dominator_tree_depth']} | {dm['pr_sized_islands_25_to_600_count']} maximal PR islands (25..600 files)"
    )
    print("  Top 10 PR-Ready Encapsulated Islands (1-Vertex-Cut, 0 interior side-entrances):")
    for isl in dm["top_20_pr_ready_islands"][:10]:
        print(
            f"    files={isl['dominated_files']:4d} (lines={isl['dominated_lines']:7,d}) | ext_deps={isl['external_foundation_deps']:3d} | "
            f"{isl['gateway_module']}"
        )
    print("  Top 5 Encapsulated Islands Strictly Inside AbelianVarietyTheory (H9):")
    for isl in dm["top_12_avt_h9_encapsulated_islands"][:5]:
        print(
            f"    files={isl['dominated_files']:4d} (lines={isl['dominated_lines']:7,d}) | ext_deps={isl['external_foundation_deps']:3d} | "
            f"{isl['gateway_module']}"
        )


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Dominator PR Target Finder & Challenge.lean 5-Theorem Audit")
    parser.add_argument("--repo", help="Path to local git clone of leanprover-community/riemann-hypothesis")
    parser.add_argument("--out", help="Write JSON snapshot to this path")
    parser.add_argument("--snapshot", help="Read existing JSON snapshot")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    if args.repo:
        data = extract_dominators_and_challenge(args.repo)
        if args.out:
            Path(args.out).write_text(json.dumps(data, indent=2), encoding="utf-8")
    elif args.snapshot:
        data = json.loads(Path(args.snapshot).read_text(encoding="utf-8"))
    else:
        default_p = Path(__file__).resolve().parent / "dominator_challenge_snapshot.json"
        data = json.loads(default_p.read_text(encoding="utf-8"))

    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print_report(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
