#!/usr/bin/env python3
"""
lean_ast_ingestor.py — Stage 1: Lean 4 Source & Export-JSON Blueprint Ingestor

Ingests Lean 4 `.lean` blueprint files (such as `rh_blueprint.lean`) AND
`lean4export`-style JSON kernel expression DAGs (`forallE` / `Exists` / `const`),
automatically extracting:

  1. Propositional Theorem DAG & `sorry` Boundary Audit:
     Identifies every `axiom` / `theorem` / `lemma`, its upstream theorem
     dependencies, and whether its proof is complete or `sorry`-bearing.

  2. Automatic Skolem Witness Poset Extraction from Quantifier Telescopes:
     Walks every `∀ (u : T)` / `∃ (v : T)` quantifier telescope across all
     declarations: whenever `∃ (v : T)` appears inside the scope of `∀ (u : T)`
     or upstream witness `w`, emits the Skolem dependency edge `u -> v` (`v`
     depends on `u`), then runs Tarjan SCC + `(H12)` scale-uniformity checks.

  3. Structured Asymptotic Exponent Term Extraction:
     Parses `@[asymptotic: <name> | exp: <linear_expr>]` annotations on theorems
     into machine-ready linear exponent vectors for the Stage 2 Tropical MILP.

Examples:
    python lean_ast_ingestor.py f:\rh\rh_blueprint.lean
    python lean_ast_ingestor.py f:\rh\rh_blueprint.lean --export-json
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import proof_dag_decompiler as dec


@dataclass
class LeanDecl:
    name: str
    decl_kind: str  # "axiom", "theorem", "lemma"
    has_sorry: bool
    theorem_deps: List[str] = field(default_factory=list)
    quantifiers: List[Tuple[str, str, str]] = field(default_factory=list)  # (q_kind: forall|exists, var, type)
    asymptotic_specs: List[Tuple[str, str]] = field(default_factory=list)  # (term_name, expr_str)


# ---------------------------------------------------------------------------
# 1. Lean 4 `.lean` Blueprint Parser
# ---------------------------------------------------------------------------

_ASYMP_RE = re.compile(r"@\[asymptotic:\s*([A-Za-z0-9_]+)\s*\|\s*exp:\s*([^\]]+)\]")
_DECL_START_RE = re.compile(r"^(axiom|theorem|lemma)\s+([A-Za-z0-9_]+)\b(.*)$")
_QUANT_RE = re.compile(r"(∀|forall|∃|exists)\s*\(\s*([A-Za-z0-9_]+)\s*:\s*([A-Za-z0-9_]+)\s*\)")
_HYP_BINDER_RE = re.compile(r"\(\s*[A-Za-z0-9_]+\s*:\s*([A-Za-z0-9_]+)\s*\)")


def parse_lean_file_text(source: str) -> List[LeanDecl]:
    """Parse a Lean 4 blueprint file into structured `LeanDecl` nodes."""
    decls: List[LeanDecl] = []
    pending_asymp: List[Tuple[str, str]] = []

    lines = source.splitlines()
    i = 0
    n = len(lines)

    while i < n:
        raw = lines[i]
        stripped = raw.strip()

        # Check structured comment annotations preceding a declaration
        for m_as in _ASYMP_RE.finditer(stripped):
            pending_asymp.append((m_as.group(1).strip(), m_as.group(2).strip()))

        if stripped.startswith("--") or not stripped or stripped.startswith(("namespace ", "end ")):
            i += 1
            continue

        m_decl = _DECL_START_RE.match(stripped)
        if not m_decl:
            i += 1
            continue

        kind, name, rest = m_decl.groups()
        block_lines = [rest]
        i += 1
        while i < n:
            nxt = lines[i].strip()
            if not nxt and (":=" in " ".join(block_lines) or kind == "axiom"):
                break
            if nxt.startswith("--") or _DECL_START_RE.match(nxt) or nxt.startswith("end "):
                break
            block_lines.append(nxt)
            i += 1

        block = " ".join(block_lines)
        has_sorry = bool(re.search(r"\bsorry\b", block))

        # Split signature before `:=` to avoid matching proof body tokens
        sig_part = block.split(":=", 1)[0]
        # Extract quantifier telescope first so we don't confuse `(x : Real)` inside ∀/∃ with theorem deps
        quants: List[Tuple[str, str, str]] = []
        for q_tok, v_name, v_type in _QUANT_RE.findall(sig_part):
            q_norm = "forall" if q_tok in ("∀", "forall") else "exists"
            quants.append((q_norm, v_name, v_type))

        # Split header binders before the outer `:` separating binders from proposition type
        # Find top-level `:` outside parentheses
        depth = 0
        colon_idx = -1
        for idx, ch in enumerate(sig_part):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth = max(0, depth - 1)
            elif ch == ":" and depth == 0:
                colon_idx = idx
                break

        binder_part = sig_part[:colon_idx] if colon_idx >= 0 else ""
        prim_types = {"Real", "Nat", "Int", "Prop", "Type", "Cert"}
        deps = [
            t_ref for t_ref in _HYP_BINDER_RE.findall(binder_part)
            if t_ref not in prim_types
        ]

        decls.append(
            LeanDecl(
                name=name,
                decl_kind=kind,
                has_sorry=has_sorry,
                theorem_deps=deps,
                quantifiers=quants,
                asymptotic_specs=list(pending_asymp),
            )
        )
        pending_asymp.clear()

    return decls


# ---------------------------------------------------------------------------
# 2. `lean4export`-Style Kernel Expression JSON Ingestor & Exporter
# ---------------------------------------------------------------------------

def export_decls_to_lean4_json(decls: Sequence[LeanDecl]) -> Dict[str, Any]:
    """
    Serialize `LeanDecl` objects into a `lean4export`-compatible JSON AST
    where propositions are nested `forallE` and `app(Exists, lam(...))` nodes.
    """
    json_decls = []
    for d in decls:
        # Build nested kernel Expr from right to left
        expr: Dict[str, Any] = {"kind": "const", "name": "True"}
        for q_kind, var, v_type in reversed(d.quantifiers):
            if q_kind == "forall":
                expr = {"kind": "forallE", "binderName": var, "binderType": v_type, "body": expr}
            else:
                expr = {
                    "kind": "app",
                    "fn": {"kind": "const", "name": "Exists"},
                    "arg": {"kind": "lam", "binderName": var, "binderType": v_type, "body": expr},
                }
        json_decls.append({
            "name": d.name,
            "declKind": d.decl_kind,
            "hasSorry": d.has_sorry,
            "usedConsts": d.theorem_deps + (["sorryAx"] if d.has_sorry else []),
            "typeExpr": expr,
            "asymptotics": [{"name": k, "expr": v} for k, v in d.asymptotic_specs],
        })
    return {"format": "lean4export-ast-v1", "declarations": json_decls}


def ingest_lean4_json(payload: Dict[str, Any]) -> List[LeanDecl]:
    """Ingest a `lean4export`-compatible JSON AST back into `LeanDecl` objects."""
    out: List[LeanDecl] = []
    for item in payload.get("declarations", []):
        quants: List[Tuple[str, str, str]] = []
        cur = item.get("typeExpr", {})
        while isinstance(cur, dict):
            k = cur.get("kind")
            if k == "forallE":
                quants.append(("forall", cur["binderName"], cur["binderType"]))
                cur = cur.get("body")
            elif k == "app" and cur.get("fn", {}).get("name") == "Exists":
                lam = cur.get("arg", {})
                quants.append(("exists", lam["binderName"], lam["binderType"]))
                cur = lam.get("body")
            else:
                break

        used = [c for c in item.get("usedConsts", []) if c != "sorryAx"]
        has_sorry = bool(item.get("hasSorry", False) or ("sorryAx" in item.get("usedConsts", [])))
        asymp = [(a["name"], a["expr"]) for a in item.get("asymptotics", [])]
        out.append(
            LeanDecl(
                name=item["name"],
                decl_kind=item["declKind"],
                has_sorry=has_sorry,
                theorem_deps=used,
                quantifiers=quants,
                asymptotic_specs=asymp,
            )
        )
    return out


# ---------------------------------------------------------------------------
# 3. Skolem Quantifier Poset & Asymptotic Term Extraction from AST
# ---------------------------------------------------------------------------

def parse_linear_exponent_expr(
    term_name: str, hypotheses: Sequence[str], expr_str: str
) -> dec.AsymptoticTerm:
    """
    Parse a linear expression like `1.0 + 2.0*alpha_n - 1.0*alpha_k` into
    `proof_dag_decompiler.AsymptoticTerm`.
    """
    c = 0.0
    coeffs = {"alpha_m": 0.0, "alpha_k": 0.0, "alpha_n": 0.0}
    # Normalize signs so we can split into additive tokens
    norm = expr_str.replace(" ", "").replace("-", "+-")
    for tok in norm.split("+"):
        if not tok:
            continue
        if "*" in tok:
            mult_s, var_s = tok.split("*", 1)
            if var_s not in coeffs:
                raise ValueError(f"Unknown exponent variable {var_s!r} in {expr_str!r}")
            coeffs[var_s] += float(mult_s)
        elif tok in coeffs:
            coeffs[tok] += 1.0
        elif tok.startswith("-") and tok[1:] in coeffs:
            coeffs[tok[1:]] -= 1.0
        else:
            c += float(tok)

    return dec.AsymptoticTerm(
        name=term_name,
        hypotheses=tuple(hypotheses),
        c=c,
        v_m=coeffs["alpha_m"],
        v_k=coeffs["alpha_k"],
        v_n=coeffs["alpha_n"],
    )


def analyze_ingested_decls(decls: Sequence[LeanDecl]) -> Dict[str, Any]:
    """
    Build the full theorem dependency graph, extract the Skolem parameter poset
    from quantifier telescopes, check `(H12)` uniformity across `sorry` boundaries,
    and compile all `@[asymptotic]` annotations into Tropical LP terms.
    """
    decl_map = {d.name: d for d in decls}

    # 1. Extract Skolem parameter dependency graph from quantifier telescopes:
    #    Within any theorem's telescope Q_1 (v_1 : T_1), ..., Q_m (v_m : T_m):
    #    each variable v_j depends on all preceding binder variables v_1, ..., v_{j-1}
    #    whenever v_j is existentially quantified (`exists`), OR when a `forall`
    #    variable represents a domain-restricted draw (e.g. Z >= x0, x in [Z, Z^(3/2)]):
    #    in telescope order, downstream binders sit strictly after upstream witnesses!
    param_kinds: Dict[str, str] = {}
    param_deps: Dict[str, set[str]] = {}

    for d in decls:
        if d.decl_kind == "axiom" and not d.quantifiers:
            param_kinds.setdefault(d.name, "Real")
            param_deps.setdefault(d.name, set())

        seen_in_telescope: List[str] = []
        for _, var, v_type in d.quantifiers:
            param_kinds.setdefault(var, v_type)
            param_deps.setdefault(var, set()).update(seen_in_telescope)
            seen_in_telescope.append(var)

    # Map canonical aliases (`x0` -> `x0_thresh`, `Z` -> `Z_window`, `x` -> `x_scale`)
    # so `verify_quantifier_poset`'s `(H12)` uniformity check applies directly!
    alias = {"x0": "x0_thresh", "Z": "Z_window", "x": "x_scale"}
    merged_kinds: Dict[str, str] = {}
    merged_deps: Dict[str, set[str]] = {}
    for k, t in param_kinds.items():
        ak = alias.get(k, k)
        merged_kinds[ak] = t
        merged_deps.setdefault(ak, set())
    for k, ds in param_deps.items():
        ak = alias.get(k, k)
        for dep in ds:
            ad = alias.get(dep, dep)
            if ad != ak:
                merged_deps[ak].add(ad)

    param_nodes = [
        dec.ParamNode(name=k, kind=merged_kinds[k], depends_on=sorted(merged_deps[k]))
        for k in merged_kinds
    ]
    poset_report = dec.verify_quantifier_poset(param_nodes)

    # 2. Theorem-level `sorry` audit
    sorry_theorems = [d.name for d in decls if d.has_sorry]
    proved_theorems = [d.name for d in decls if not d.has_sorry and d.decl_kind != "axiom"]

    # 3. Extract asymptotic terms tagged with their theorem & transitive upstream hypotheses
    asymp_terms: List[dec.AsymptoticTerm] = []
    for d in decls:
        for t_name, expr_str in d.asymptotic_specs:
            hyps = [d.name] + [dep for dep in d.theorem_deps if dep in decl_map]
            asymp_terms.append(parse_linear_exponent_expr(t_name, hyps, expr_str))

    return {
        "total_declarations": len(decls),
        "sorry_bearing_declarations": sorry_theorems,
        "proved_declarations": proved_theorems,
        "quantifier_poset": poset_report,
        "extracted_asymptotic_terms": [
            {
                "name": t.name,
                "hypotheses": list(t.hypotheses),
                "c": t.c,
                "v_m": t.v_m,
                "v_k": t.v_k,
                "v_n": t.v_n,
            }
            for t in asymp_terms
        ],
        "_raw_terms": asymp_terms,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 1: Lean 4 Blueprint & Export-JSON Ingestor")
    parser.add_argument("path", nargs="?", default=str(Path(__file__).resolve().parent / "rh_blueprint.lean"))
    parser.add_argument("--export-json", action="store_true", help="Emit lean4export-compatible JSON AST")
    args = parser.parse_args(argv)

    text = Path(args.path).read_text(encoding="utf-8")
    decls = parse_lean_file_text(text)
    if args.export_json:
        print(json.dumps(export_decls_to_lean4_json(decls), indent=2))
    else:
        rep = analyze_ingested_decls(decls)
        rep.pop("_raw_terms", None)
        print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
