# Auditing the Alpoge–Claude Riemann Hypothesis Proof: The Local-to-Global Step

**Explicit derivation, structural resolution of §6.3, `(H0)–(H12)` Lean 4 Rosetta Stone, Dominator PR Finder, and 76-test + 11-theorem Lean 4 verification suite.**

[![Verify Mathematical & Computational Suite](https://github.com/seedhunterai/riemann-hypothesis-local-to-global/actions/workflows/verify.yml/badge.svg)](https://github.com/seedhunterai/riemann-hypothesis-local-to-global/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Companion monograph and reproducible exact-arithmetic codebase examining the local-to-global mechanism in the public reconstruction of the Alpoge–Claude proof of the Riemann hypothesis, resolving the structural questions of §6.3, isolating the non-abelian $`GL_2`$ barrier, and analyzing the polyhedral structure of window-free Quasi-RH reductions.

- **Full Monograph:** [`LOCAL_TO_GLOBAL.md`](LOCAL_TO_GLOBAL.md)
- **Live Interactive Explorer:** **[https://seedhunterai.github.io/riemann-hypothesis-local-to-global/](https://seedhunterai.github.io/riemann-hypothesis-local-to-global/)** *(source: [`docs/index.html`](docs/index.html))*
- **Zero-`sorry` Lean 4 Proofs:** [`lean/LocalToGlobal/Theorems.lean`](lean/LocalToGlobal/Theorems.lean) *(11 kernel-checked theorems including Conductor-Aspect Quasi-GRH, zero external dependencies, `lake build` in `~1s`)*
  - **Part I (§§0–8):** Explicit derivation of the zero-side Mellin/Landau window reduction (`[D]`), the Stuhler–Grayson slope inequality engine (`[D]`), the main inequality `(M)` (`[Q]`), and the audit of `(H0)–(H12)` (`[H]`).
  - **Part II (§§9–13):** Proved theorems (`[D]`) and exact arbitrary-precision computational certificates (`[C]`) on:
    1. **§9 — Determinant Identities vs. Padé Noise (§6.3.2):** Exact triangular basis invariance (Theorem 9.1 `[D]`), $`208`$-digit unshifted hyperfactorial identity, prime-certified split-prime counterexamples refuting single-determinant local factorization at shifted scales (Proposition 9.2 `[D]`), $`51`$-scale Padé annihilation via the zeroth Fitting ideal (`[C]`), Cauchy–Binet kernel trap & Malmsten–Kronecker lemniscate identity (Proposition 9.3 `[D]`), and rank-deficient ($`r < D`$) modular Hermite Normal Form spans (`[C]`).
    2. **§10 — Non-Abelian $`GL_2`$ Barriers:** Honda integrality, ordinary-prime slope collapse, $`100\%`$ purely $`p`$-adic Frobenius trace recovery mod $`p^2`$ without $`\mathbb{F}_p`$ point counting, and $`a_p`$-free $`2\times 2`$ Hecke–Wronskians (Theorem 10.1 `[D]`).
    3. **§11 — Projective Towers (§6.3.1):** Exact horizontal pullback and $`2`$-adic Iwasawa distribution pushforward $`\mu_{2N}(c_+)+\mu_{2N}(c_-)=\mu_N(c)`$ across dyadic levels vs. vertical infinite linear rank over $`\mathbb{Q}`$ (Theorem 11.1 `[D]`).
    4. **§12 — Polyhedral Decompilation & Window-Free Quasi-RH:** Skolem quantifier-poset verification, hidden-constraint LP gap, and exact rational Farkas certification of $`\theta^*(\tau_k)=\max(\tfrac12,\,2-\tfrac{\tau_k}{3})`$ on the blueprint model (Proposition 12.1 `[D]`), plus the Guth–Maynard (`30/13`) zero-density polyhedron saddle (`[C]`).
    5. **Live Upstream Lean 4 DAG Audit (`upstream_lean_auditor.py`):** Direct streaming audit of all **58,880 `.lean` files** (`13.93M` lines, `931 MB`) in [`leanprover-community/riemann-hypothesis`](https://github.com/leanprover-community/riemann-hypothesis) (`ce1560fb`), verifying the exact acyclic import DAG from `src/Main.lean` (`58,858` reachable modules, `359,031` edges, max chain depth `95`, `0` reachable `sorry`s, `0` custom axioms) and mapping the codebase volume onto `(H0)–(H12)`.
    6. **Verified 16-Row Rosetta Stone & `H9` Boundary Cut-Set (`upstream_cutset_rosetta.py`):** Exact file/line/signature permalinks mapping `Main`, `Lemmas 1.1–1.2`, `Lemma 2.1`, `(H0)–(H12)`, and `Estimates B & C` onto `ce1560fb`, plus the exact cross-boundary vertex cut-set of `Development.AbelianVarietyTheory` (`10,292` direct gateway modules = **`39.95%`** of its `25,765` modules across `33,148` incoming edges).
    7. **Dominator-Tree `"PR Target Finder"` & `Challenge.lean` Conductor-Aspect Quasi-GRH (`upstream_dominator_and_challenge.py`):** Exact DAG dominator decomposition of all `58,858` reachable modules identifying **`4,413` certified `0`-side-entrance refactoring islands** (`155` PR-sized cones of `25–600` files, including `1,268`-file encapsulated islands inside `H9`), Nov 20 $\to$ `ce1560fb` clumping/debt telemetry (`175,187` $\to$ `58,880` files, `627` `set_option backward` files), the `src/Challenge.lean` 5-theorem inversion (`riemannHypothesis := Claim.zeta LFunction_ne_zero_of_half_lt_re`), and the Lean 4 kernel-verified Conductor-Aspect ($`q = x^{\alpha_q}`$) Quasi-GRH closed form $`\theta^*(\tau_k,\alpha_q)=\max\!\left(\tfrac{1+\alpha_q}{2},\,2+\alpha_q-\tfrac{\tau_k}{3}\right)`$.

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

> **Key takeaway:** The local Lubin–Tate prime-divisibility engine (`H6–H7`) accounts for **under 0.8%** of the live Lean proof (`111k` lines), while the global abelian-variety/cohomology denominator machinery (`H9`) alone accounts for **46.16%** (`6.46M` lines). Full JSON snapshot in [`upstream_audit_snapshot.json`](upstream_audit_snapshot.json).

### Verified 16-Row `(H0)–(H12)` Rosetta Stone & `H9` Cut-Set (`ce1560fb`)

> **The `H9` Mycelial Mesh Discovery (`rosetta_cutset_snapshot.json`):** Unlike human-authored libraries with a narrow API facade ($`<1\%`$ of modules), out of **25,765** reachable modules inside `Development.AbelianVarietyTheory` (`6,143,414` lines), **10,292 distinct modules (`39.95%`)** are directly imported from outside across **33,148** cross-boundary edges by **5,653** external callers (led by `MatchingOfParts` [`247`], `WeilTensorModels` [`225`], and `MkMemAdmissibleAtIffCocycle` [`192`]). Refactoring or replacing `(H9)` therefore cannot proceed by folder-level import cutting; it must sever at the single downstream theorem interface `UpperTheory.One.vwin_of_row` ([`VwinOfRow.lean#L1230`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/AutomorphicVarieties/ModularityArgument/ThresholdArguments/DenseThresholdCases/VwinOfRow.lean#L1230)).

| Monograph Item | Section | Verified Lean 4 Permalink (`@ ce1560fb`) | Declaration Signature Excerpt | Tropical LP / Role |
|---|---|---|---|---|
| **Main Headline** | §0 | [`Main.lean#L12`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Main.lean#L12) | `theorem riemannHypothesis : RiemannHypothesis` | `ROOT_GOAL` |
| **Lemma 1.1 & 1.2** | §1.2–1.3 | [`ZeroFreeAllOfTopG.lean#L1151`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/CharacterBoxes/ZeroFreeAllOfTopG.lean#L1151) | `theorem zeroFree_of_rpow_window` | `ANALYTIC_ZERO_SIDE` |
| **Lemma 2.1** | §2.3 | [`FrameHeightOfDualSlope.lean#L461`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/LocalToGlobal/ExponentGrids/ExponentGridsProofs/FrameHeightOfDualSlope.lean#L461) | `theorem grid_function_of_lemma_A` | `LINEAR_ALGEBRA_CORE` |
| **`(H0)`** | §3.1, §5.2 | [`DataParams.lean#L53`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/LocalRingCardinality/DataParams.lean#L53) | `structure ParameterChoice' (D : D7) (L : SizeData₂ D)` | `CRITICAL_IN_QUASI_RH` ($`\lambda^*=2/3`$) |
| **`(H1)`** | §3.2(b) | [`EveryCharacterOfLiteratureOfVWin.lean#L958`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean#L958) | `published_results_analytic` | `SLACK_AT_4/5 (+0.30)` |
| **`(H2)`** | §6.1 | [`EveryCharacterOfLiteratureOfVWin.lean#L971`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean#L971) | `published_results_local : GlobalTheory.NLink.Hypotheses₂` | `SLACK_AT_4/5 (+0.30)` |
| **`(H3)`** | §6.1 | [`DataNLoc.lean#L1214`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/SharpRestriction/LevelLetterCocycles/DataNLoc.lean#L1214) | `structure LocHW (l : ℕ) [Fact l.Prime] (P : LocPrime l)` | `CRITICAL_IN_QUASI_RH` |
| **`(H4)`** | §3.3, §6.1 | [`GenusMassVolumes.lean#L365`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/SharpRestriction/GlobalSchurBounds/GeneralLinearActions/GenusMassVolumes.lean#L365) | `def SdHerm ... (ψR : R → ℂ) (n : ℕ)` | `SLACK_AT_4/5 (+0.30)` |
| **`(H5)`** | §3.2(c) | [`EveryCharacterOfLiteratureOfVWin.lean#L1065`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean#L1065) | `published_results_arithmetic : ...PublishedHodgeBundleTheorems` | `SLACK_AT_4/5 (+0.30)` |
| **`(H6)`** | §3.2(a), §9 | [`AbsNormEqRelIndexMul.lean#L21`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/LocalToGlobal/CMDeterminantIdeals/CMIdealValuations/AbsNormEqRelIndexMul.lean#L21) | `theorem absNorm_eq_relIndex_mul` | `SLACK_AT_4/5 (+0.30)` |
| **`(H7)`** | §3.2(a) | [`FinalLL.lean#L20`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/ThetaLogCoordinates/BihomogeneousPairs/FinalLL.lean#L20) | `theorem levelLattice : Formalism.GM.M.LevelLattice` | `LOCAL_TO_GLOBAL_BRIDGE` |
| **`(H8)`** | §2.4 | [`DefsPl.lean#L18`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/HarderNarasimhan/DefsPl.lean#L18) | `def PlaceSplit : Prop` | `LINEAR_ALGEBRA_FILTRATION` |
| **`(H9)`** | §3.2(c), §5.3 | [`VwinOfRow.lean#L1230`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/AutomorphicVarieties/ModularityArgument/ThresholdArguments/DenseThresholdCases/VwinOfRow.lean#L1230) | `theorem UpperTheory.One.vwin_of_row ... : UpperTheory.Top.VWin` | `WINDOW_GENERATOR (46.2%)` |
| **`(H10) & (H11)`** | §2.4, §12 | [`Assemble.lean#L28`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Formalism/LocalToGlobal/PairingIndexBounds/HalvesRecursion/BoundedHalfModels/Assemble.lean#L28) | `theorem halvesRec_of_facts` | `CRITICAL_IN_QUASI_RH` ($`\lambda^*=1/3`$) |
| **Estimates B & C** | §4 | [`EveryCharacterOfLiteratureOfVWin.lean#L1083`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean#L1083) | `every_character_of_literature_of_vWin` | `AUX_PRIME_SUPPLY` |
| **`(H12)`** | §5.1 | [`ZeroFreeAllOfTopG.lean#L1259`](https://github.com/leanprover-community/riemann-hypothesis/blob/ce1560fba43e59629abdf3c27a279218fad9f77a/src/Development/AutomorphicArithmetic/CharacterBoxes/ZeroFreeAllOfTopG.lean#L1259) | `theorem zeroFreeAll_of_topG ... : ZeroFreeAll` | `SCALE_UNIFORMITY_QUANTIFIER` |

### Contributor `"PR Target Finder"` (Dominator Islands) & `Challenge.lean` Inversion

Scanning the exact immediate-dominator tree $`T_{\mathrm{dom}}`$ from `src/Main.lean` ([`dominator_challenge_snapshot.json`](dominator_challenge_snapshot.json)) isolates **`4,413` certified `1`-vertex-cut islands** (`0` interior side-entrances across all `359,031` edges), plus tracks the community cleanup from the Nov 20 baseline (`175,187` files / `17.37M` lines $\to$ `58,880` files / `13.93M` lines = **`-66.39%` files**, via `12,024` auto-clumped files holding `63,816` clump members; `627` `set_option backward` and `3` `exact_by_kernel` debt files remaining):

- **Macro Dominator Spine & `Challenge.lean` Inversion:** In `src/Main.lean`, `riemannHypothesis := Claim.zeta LFunction_ne_zero_of_half_lt_re` — `MainConclusion` proves `DirichletLNonvanishingRightOfHalf` across **`58,857` files (`99.998%`)** first, and $`\zeta(s)`$ specializes $`N=1, \chi=1`$ in `ClaimStrip.lean` (`4`-file cone). Over half the entire repo (**`30,024` files, `51.01%`**) hangs strictly beneath `EveryCharacterOfLiteratureOfVWin`, and **`15,677` files (`26.64%`)** hang strictly beneath `SiegelEisensteinSecondTermHypotheses`.
- **Sample Atomic PR Cones (`side_entrances = 0`):** `ToroidalCompactificationFanLogKodairaSpencerIso` (**`1,268` files / `353k` lines** strictly inside `H9`!), `MainTheoremPublishedInputsThird` (`587` files), `SphAtTw` (`580` files, only `229` ext. deps), `PoincareLelongInv` (`545` files proving `poincareLelong`, `267` ext. deps), `EtSideGood` (`447` files, `186` ext. deps), and `MoriwakiBigTheorem.Outright` (`272` files proving `moriwakiBig_outright`).
- **Conductor-Aspect ($`q = x^{\alpha_q}`$) Window-Free Quasi-GRH (`Theorems.lean` §5 `[D]`):** Extending `QuasiRHModel` to primitive Dirichlet $`L(s,\chi)`$ of conductor $`q = x^{\alpha_q}`$ yields the exact rational Farkas closed form $`\theta^*(\tau_k,\alpha_q) = \max\!\left(\tfrac{1+\alpha_q}{2},\, 2 + \alpha_q - \tfrac{\tau_k}{3}\right)`$, hitting the tropical convexity barrier $`\theta^* \ge 1`$ at $`\alpha_q = 1`$ ($`q = x`$).

---

## Epistemic Discipline

Every claim in [`LOCAL_TO_GLOBAL.md`](LOCAL_TO_GLOBAL.md) carries an explicit status tag:
- **`[D]`** — Proved unconditionally in the text for all parameter values.
- **`[C]`** — Certified by exact integer or exact `Fraction` arithmetic in this repository (`76` automated tests). Finite-range computational findings are never overstated as universal theorems without proof.
- **`[Q]` / `[H]`** — Quoted statements and unproved black-box hypotheses from M. Rubinstein's public reconstruction (*working draft v2*, 22 November 2026). This repository does **not** prove `(H0)–(H12)` and does **not** claim an independent proof of RH.

---

## Quickstart

```bash
pip install -r requirements.txt
pytest -v
```

Running `pytest -v` executes all **76 verification tests** across 13 test modules in ~45 seconds.

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

# 16-row (H0)-(H12) Rosetta Stone & AbelianVarietyTheory boundary cut-set audit
python upstream_cutset_rosetta.py --snapshot rosetta_cutset_snapshot.json
python upstream_cutset_rosetta.py --repo /path/to/riemann-hypothesis --out rosetta_cutset_snapshot.json

# Dominator-tree PR Target Finder, Challenge.lean cone split, & Conductor Quasi-GRH
python upstream_dominator_and_challenge.py --snapshot dominator_challenge_snapshot.json
python upstream_dominator_and_challenge.py --repo /path/to/riemann-hypothesis --out dominator_challenge_snapshot.json
```
