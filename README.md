# The Local-to-Global Step: Monograph & Computational Verification Suite

[![Verify Mathematical & Computational Suite](https://github.com/seedhunterai/rh-local-to-global/actions/workflows/verify.yml/badge.svg)](https://github.com/seedhunterai/rh-local-to-global/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Companion monograph and reproducible exact-arithmetic codebase examining the local-to-global mechanism in the public reconstruction of the Alpoge–Claude proof of the Riemann hypothesis, resolving the structural questions of §6.3, isolating the non-abelian $GL_2$ barrier, and analyzing the polyhedral structure of window-free Quasi-RH reductions.

- **Full Monograph:** [`LOCAL_TO_GLOBAL.md`](LOCAL_TO_GLOBAL.md)
  - **Part I (§§0–8):** Explicit derivation of the zero-side Mellin/Landau window reduction (`[D]`), the Stuhler–Grayson slope inequality engine (`[D]`), the main inequality `(M)` (`[Q]`), and the audit of `(H0)–(H12)` (`[H]`).
  - **Part II (§§9–13):** Proved theorems (`[D]`) and exact arbitrary-precision computational certificates (`[C]`) on:
    1. **§9 — Determinant Identities vs. Padé Noise (§6.3.2):** Exact triangular basis invariance (Theorem 9.1 `[D]`), $208$-digit unshifted hyperfactorial identity, prime-certified split-prime counterexamples refuting single-determinant local factorization at shifted scales (Proposition 9.2 `[D]`), $51$-scale Padé annihilation via the zeroth Fitting ideal (`[C]`), Cauchy–Binet kernel trap & Malmsten–Kronecker lemniscate identity (Proposition 9.3 `[D]`), and rank-deficient ($r < D$) modular Hermite Normal Form spans (`[C]`).
    2. **§10 — Non-Abelian $GL_2$ Barriers:** Honda integrality, ordinary-prime slope collapse, $100\%$ purely $p$-adic Frobenius trace recovery mod $p^2$ without $\mathbb{F}_p$ point counting, and $a_p$-free $2\times 2$ Hecke–Wronskians (Theorem 10.1 `[D]`).
    3. **§11 — Projective Towers (§6.3.1):** Exact horizontal pullback and $2$-adic Iwasawa distribution pushforward $\mu_{2N}(c_+)+\mu_{2N}(c_-)=\mu_N(c)$ across dyadic levels vs. vertical infinite linear rank over $\mathbb{Q}$ (Theorem 11.1 `[D]`).
    4. **§12 — Polyhedral Decompilation & Window-Free Quasi-RH:** Skolem quantifier-poset verification, hidden-constraint LP gap, and exact rational Farkas certification of $\theta^*(\tau_k)=\max(\tfrac12,\,2-\tfrac{\tau_k}{3})$ on the blueprint model (Proposition 12.1 `[D]`), plus the Guth–Maynard (`30/13`) zero-density polyhedron saddle (`[C]`).

---

## Epistemic Discipline

Every claim in [`LOCAL_TO_GLOBAL.md`](LOCAL_TO_GLOBAL.md) carries an explicit status tag:
- **`[D]`** — Proved unconditionally in the text for all parameter values.
- **`[C]`** — Certified by exact integer or exact `Fraction` arithmetic in this repository (`60` automated tests). Finite-range computational findings are never overstated as universal theorems without proof.
- **`[Q]` / `[H]`** — Quoted statements and unproved black-box hypotheses from M. Rubinstein's public reconstruction (*working draft v2*, 22 November 2026). This repository does **not** prove `(H0)–(H12)` and does **not** claim an independent proof of RH.

---

## Quickstart

```bash
pip install -r requirements.txt
pytest -v
```

Running `pytest -v` executes all **60 verification tests** across 10 test modules in ~45 seconds.

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
```
