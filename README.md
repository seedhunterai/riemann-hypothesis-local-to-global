# Auditing the Alpoge–Claude Riemann Hypothesis Proof: The Local-to-Global Step

**Explicit derivation, structural resolution of §6.3, and 60-test verification suite for the Riemann zeta function and Lean blueprint.**

[![Verify Mathematical & Computational Suite](https://github.com/seedhunterai/riemann-hypothesis-local-to-global/actions/workflows/verify.yml/badge.svg)](https://github.com/seedhunterai/riemann-hypothesis-local-to-global/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Companion monograph and reproducible exact-arithmetic codebase examining the local-to-global mechanism in the public reconstruction of the Alpoge–Claude proof of the Riemann hypothesis, resolving the structural questions of §6.3, isolating the non-abelian $`GL_2`$ barrier, and analyzing the polyhedral structure of window-free Quasi-RH reductions.

- **Full Monograph:** [`LOCAL_TO_GLOBAL.md`](LOCAL_TO_GLOBAL.md)
- **Live Interactive Explorer:** **[https://seedhunterai.github.io/riemann-hypothesis-local-to-global/](https://seedhunterai.github.io/riemann-hypothesis-local-to-global/)** *(source: [`docs/index.html`](docs/index.html))*
  - **Part I (§§0–8):** Explicit derivation of the zero-side Mellin/Landau window reduction (`[D]`), the Stuhler–Grayson slope inequality engine (`[D]`), the main inequality `(M)` (`[Q]`), and the audit of `(H0)–(H12)` (`[H]`).
  - **Part II (§§9–13):** Proved theorems (`[D]`) and exact arbitrary-precision computational certificates (`[C]`) on:
    1. **§9 — Determinant Identities vs. Padé Noise (§6.3.2):** Exact triangular basis invariance (Theorem 9.1 `[D]`), $`208`$-digit unshifted hyperfactorial identity, prime-certified split-prime counterexamples refuting single-determinant local factorization at shifted scales (Proposition 9.2 `[D]`), $`51`$-scale Padé annihilation via the zeroth Fitting ideal (`[C]`), Cauchy–Binet kernel trap & Malmsten–Kronecker lemniscate identity (Proposition 9.3 `[D]`), and rank-deficient ($`r < D`$) modular Hermite Normal Form spans (`[C]`).
    2. **§10 — Non-Abelian $`GL_2`$ Barriers:** Honda integrality, ordinary-prime slope collapse, $`100\%`$ purely $`p`$-adic Frobenius trace recovery mod $`p^2`$ without $`\mathbb{F}_p`$ point counting, and $`a_p`$-free $`2\times 2`$ Hecke–Wronskians (Theorem 10.1 `[D]`).
    3. **§11 — Projective Towers (§6.3.1):** Exact horizontal pullback and $`2`$-adic Iwasawa distribution pushforward $`\mu_{2N}(c_+)+\mu_{2N}(c_-)=\mu_N(c)`$ across dyadic levels vs. vertical infinite linear rank over $`\mathbb{Q}`$ (Theorem 11.1 `[D]`).
    4. **§12 — Polyhedral Decompilation & Window-Free Quasi-RH:** Skolem quantifier-poset verification, hidden-constraint LP gap, and exact rational Farkas certification of $`\theta^*(\tau_k)=\max(\tfrac12,\,2-\tfrac{\tau_k}{3})`$ on the blueprint model (Proposition 12.1 `[D]`), plus the Guth–Maynard (`30/13`) zero-density polyhedron saddle (`[C]`).
    5. **Live Upstream Lean 4 DAG Audit (`upstream_lean_auditor.py`):** Direct streaming audit of all **58,880 `.lean` files** (`13.93M` lines, `931 MB`) in [`leanprover-community/riemann-hypothesis`](https://github.com/leanprover-community/riemann-hypothesis) (`ce1560fb`), verifying the exact acyclic import DAG from `src/Main.lean` (`58,858` reachable modules, `359,031` edges, max chain depth `95`, `0` reachable `sorry`s, `0` custom axioms) and mapping the codebase volume onto `(H0)–(H12)`.

---

## Live Upstream Audit (`leanprover-community/riemann-hypothesis @ ce1560fb`)

| Architectural Role (`(H0)–(H12)`) | `.lean` Files | Share of `src/` | Reachable from `Main` | Total Lines |
|---|---:|---:|---:|---:|
| **`H9` — Abelian Varieties & Galois/Étale/Profinite Cohomology** | **27,174** | **46.16%** | 27,172 | **6,457,976** |
| **`H1, H4, H5` — Heights, Theta Lifts, Periods & Mass Formulas** | 9,401 | 15.97% | 9,401 | 1,728,237 |
| **`Est. B & C` — Automorphic Spectrum, Weil Bound & $`GL_3`$ Mean Squares** | 7,537 | 12.80% | 7,533 | 1,915,491 |
| **`H12` — Main Reduction, Character Boxes & Windows (`vWin`)** | 7,068 | 12.01% | 7,063 | 1,984,819 |
| **`H8, H10` — Local-to-Global, Projective Lattices & Dual Slopes** | 5,428 | 9.22% | 5,427 | 1,356,228 |
| **Other Supporting Modules** | 1,561 | 2.65% | 1,560 | 313,098 |
| **`H6, H7` — Local Formal Groups & Lubin–Tate Divisibility** | **415** | **0.70%** | 415 | **111,939** |
| **`H0` — Global Fields & Foundations** | 287 | 0.49% | 287 | 65,326 |
| **Total (`src/`)** | **58,871** | **100.00%** | **58,858** | **13,933,114** |

> **Key takeaway:** The local Lubin–Tate prime-divisibility engine (`H6–H7`) accounts for **under 0.8%** of the live Lean proof (`111k` lines), while the global abelian-variety/cohomology denominator machinery (`H9`) alone accounts for **46.35%** (`6.46M` lines). Full JSON snapshot in [`upstream_audit_snapshot.json`](upstream_audit_snapshot.json).

---

## Epistemic Discipline

Every claim in [`LOCAL_TO_GLOBAL.md`](LOCAL_TO_GLOBAL.md) carries an explicit status tag:
- **`[D]`** — Proved unconditionally in the text for all parameter values.
- **`[C]`** — Certified by exact integer or exact `Fraction` arithmetic in this repository (`64` automated tests). Finite-range computational findings are never overstated as universal theorems without proof.
- **`[Q]` / `[H]`** — Quoted statements and unproved black-box hypotheses from M. Rubinstein's public reconstruction (*working draft v2*, 22 November 2026). This repository does **not** prove `(H0)–(H12)` and does **not** claim an independent proof of RH.

---

## Quickstart

```bash
pip install -r requirements.txt
pytest -v
```

Running `pytest -v` executes all **64 verification tests** across 11 test modules in ~45 seconds.

### Key Standalone Reports

```bash
# Baseline & 4 use-case diagnostics (Rubinstein's small table, double-zero blindness, Sato-Tate)
python rh_audit.py report

# GL_2 local-to-global analytic rank & Weil-Mestre GRH ceiling solver (ranks 0, 1, 2, 3)
python gl2_local_to_global.py --benchmark

# §9: Unshifted 208-digit hyperfactorial vs. shifted split-prime Padé obstruction
python cm_determinant_lab.py report

# §9.3 & §10.3: 51-scale Fitting-ideal Padé annihilation (324 -> 201 digits) & GL_2 Wronskians
python fitting_ideal_lab.py report

# §9.3 & §5.3: Gauss-Kronecker 10^-19 check, Bezout law, & ell-adic twist ball scanner
python frontier6_bezout_and_centres.py report

# §2.4 & §9.3: Rank-deficient (r = D - 1 < D) multi-word modular HNF simulator
python frontier7_hnf_slope_engine.py report

# §10: Honda integrality, ordinary slope collapse, & mod-p^2 ASD trace recovery
python gl2_divisibility_lab.py report

# §11: Dyadic Iwasawa measure exactness (N = 8..256) vs. vertical infinite rank
python projective_tower_lab.py report

# §12: Lean blueprint AST ingestion, Quasi-RH Q-Farkas dichotomy, & Guth-Maynard 30/13 saddle
python quasi_rh_extractor.py report

# Live audit of leanprover-community/riemann-hypothesis (snapshot or fresh local git clone)
python upstream_lean_auditor.py --snapshot upstream_audit_snapshot.json
python upstream_lean_auditor.py --repo /path/to/riemann-hypothesis --out upstream_audit_snapshot.json
```
