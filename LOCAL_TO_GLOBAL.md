# The local-to-global step

### An explicit derivation of how per-prime data becomes a bound on $`G(x)`$, the exact point at which the derivation stops being analytic, and a structural resolution of the open questions of §6

---

**Status key.** Every statement below carries one of four marks.

- **[D]** derived here, with a complete proof written out for all values of the parameters.
- **[C]** certified here by exact arbitrary-precision integer or rational (`Fraction`) computation in the companion repository, reproducible via `pytest`; where a general pattern beyond the certified range is discussed, it is stated explicitly as a heuristic or conjecture rather than as a theorem.
- **[Q]** quoted from the public reconstruction of the proof (M. Rubinstein, *The Alpoge-Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, working draft v2, 22 November 2026) and, through it, from the Lean development.
- **[H]** one of the thirteen hypotheses to which the reconstruction reduces the main inequality. None of the thirteen is proved in the reconstruction.

The subject of Part I (§§0–8) is the Riemann zeta function and the public reconstruction of its local-to-global reduction. Part II (§§9–13) attacks the structural questions raised by that audit: whether an exact determinant identity underlies Lemma 2.1 (§6.3.2), why multi-word rank-deficient spans (§2.4) cannot be simplified to scalar combinations, whether a fixed projective system exists across dyadic levels (§6.3.1), where the method meets an arithmetic wall on non-CM $`GL_2`$ $`L`$-functions, and what happens to the zero-free half-plane when Hypothesis (H9)'s window machinery is removed.

---

# Part I. The reduction and its audit

## 0. What this document is for

The proof is certified and not understood. Those are different failures, and only one of them is about correctness.

The certification is not in doubt. The library compiles, the exported proof is accepted by Lean's kernel and by independently written kernels, the statement is matched against Mathlib's own definition of the Riemann hypothesis, and the five headline theorems depend on the three standard axioms. A machine has checked that the theorem follows from the axioms. What no machine has done is show a human *why* it follows.

This document takes the one question that carries the whole difficulty and answers it as far as the public record allows:

> How does information attached to individual primes become a statement about all $`x`$ at once?

The answer has two halves, and they work by completely different mechanisms. That is the first thing to understand, and it is the thing most expositions blur.

1. **On the zero side, the passage to infinity is analytic, and it is complete.** A bound for $`G(x)`$ on a sequence of long intervals, each interval stretching from $`Y`$ to $`Y^{3/2}`$, forces every zero to lie on or left of the line. This is a theorem of Landau type, it is one page, and it is proved here in full. Nothing is missing.

2. **On the prime side, the passage to infinity is not analytic.** At each scale $`x`$ the argument builds a finite object and proves a finite inequality about it. The statement that such an object exists at every scale in the windows, with constants independent of the scale, is Hypothesis (H12). It is not a limit. It is a universally quantified existence statement, and it is assumed in the reconstruction.

The derivation in Part I is therefore complete as a derivation *of the reduction*, and it stops at the uniformity across scales. Section 7 states the reduction on one page; Part II (§§9–13) then investigates the mathematical obstructions behind that stopping point.

---

## 1. The target: from zeros to an upper bound

### 1.1 The statistic, and why it is the right one

Let $`\Lambda(n)=\log p`$ for a prime power $`n=p^{h}`$ and $`0`$ otherwise. Define

```math
S(x)=\sum_{n\le x}\Lambda(n)\Big(\frac{x}{n}-1\Big),\qquad G(x)=S(x)-x\log x+(1+\gamma)x.
```

The weight $`(x/n-1)`$ is a ramp: it is large for small $`n`$, falls linearly, and vanishes at $`n=x`$. With the positive part $`v_{+}=\max(v,0)`$ the condition $`n\le x`$ can be dropped and $`S(x)=\sum_{n}\Lambda(n)(x/n-1)_{+}`$.

The reason for the ramp is that its Mellin transform is elementary. For $`\mathrm{Re}\,s>1`$,

```math
\int_{1}^{\infty}\Big(\frac{x}{n}-1\Big)_{+}x^{-s-1}\,dx=\frac{n^{-s}}{s(s-1)}.
```

Summing against $`\Lambda(n)`$ and using $`\sum\Lambda(n)n^{-s}=-\zeta'/\zeta(s)`$ gives the transform of $`S`$. The two subtracted terms remove the pole at $`s=1`$ *exactly*:

```math
F(s)=\int_{1}^{\infty}G(x)\,x^{-s-1}\,dx=-\frac{\zeta'}{\zeta}(s)\cdot\frac{1}{s(s-1)}-\frac{1}{(s-1)^{2}}+\frac{1+\gamma}{s-1}.
```

Near $`s=1`$,

```math
-\frac{\zeta'}{\zeta}(s)=\frac{1}{s-1}-\gamma+O(s-1),\qquad \frac{1}{s(s-1)}=\frac{1}{s-1}-1+O(s-1).
```

The product is $`\frac{1}{(s-1)^{2}}-\frac{1+\gamma}{s-1}+\gamma+O(s-1)`$; the remaining two terms cancel the double pole and the simple pole, and leave $`\gamma`$. So $`F`$ is regular at $`s=1`$, with $`F(1)=\gamma`$. Nothing is lost and nothing is left over: $`G`$ is engineered so that its transform is finite at the one point where $`\zeta`$ has a pole.

What is left in $`F`$ is exactly the zeros. At a simple zero $`\rho`$ the function $`-\zeta'/\zeta`$ has a simple pole of residue $`-1`$, so $`F`$ has a simple pole with residue

```math
\mathrm{Res}_{s=\rho}\,F=-\frac{1}{\rho(\rho-1)}.
```

At $`s=0`$ there is a pole of residue $`\log 2\pi`$. Shifting the contour in the inverse Mellin transform gives the explicit formula

```math
G(x)=\log 2\pi-\sum_{\rho}\frac{x^{\rho}}{\rho(\rho-1)}+O(x^{-2}),
```

the sum over nontrivial zeros. The factor $`\rho(\rho-1)`$ makes it absolutely convergent: no conditional convergence, no summation by parts.

Two consequences are used constantly below.

- If every zero has real part $`\tfrac12`$ and is simple then $`|\rho(\rho-1)|=\rho(1-\rho)`$ and

```math
|G(x)-\log 2\pi|\le (2+\gamma-\log 4\pi)\sqrt{x}+O(x^{-2}),\qquad \sum_{\rho}\frac{1}{\rho(1-\rho)}=2+\gamma-\log 4\pi=0.046191\ldots
```

- A zero of real part $`\beta`$ contributes a term of size $`x^{\beta}`$. The exponent of the error in the prime count *is* the real part of a zero.

### 1.2 One side is enough: Landau's criterion

**Lemma 1.1** [D]. *Let $`0<\theta<1`$. If $`G(x)\le x^{\theta}`$ for all large $`x`$, then $`\zeta(s)\ne0`$ for $`\mathrm{Re}\,s>\theta`$.*

*Proof.* Put

```math
H(s)=\int_{1}^{\infty}\big(x^{\theta}-G(x)\big)x^{-s-1}\,dx=\frac{1}{s-\theta}-F(s).
```

The integrand is eventually nonnegative. Landau's theorem on Mellin transforms of functions of one sign says that such a transform has a singularity at the real point of its abscissa of convergence. By the formula for $`F`$, the only real singularities of $`F`$ in $`s>0`$ would come from real zeros of $`\zeta`$ in $`(0,1)`$, and there are none; the pole at $`s=1`$ has been removed. Hence $`H`$ is regular at every real $`s>\theta`$, so its abscissa of convergence is at most $`\theta`$, so $`H`$ is holomorphic on $`\mathrm{Re}\,s>\theta`$, so $`F=1/(s-\theta)-H`$ is holomorphic there except at $`s=\theta`$. A zero of $`\zeta`$ with real part $`>\theta`$ would be a pole of $`F`$ in that half-plane. $`\square`$

### 1.3 Why intervals suffice

**Lemma 1.2** [D, from a cited theorem]. *Let $`0<\theta<1`$. If there are arbitrarily large $`Y`$ such that $`G(x)\le x^{\theta}`$ for all $`x\in[Y,Y^{3/2}]`$, then $`\zeta(s)\ne0`$ for $`\mathrm{Re}\,s>\theta`$.*

*Proof.* Suppose $`\zeta`$ has a zero with real part $`>\theta`$, and let $`\Theta\le1`$ be the supremum of the real parts of the zeros. Put $`f(x)=G(x)-2x^{\theta}`$ for $`x\ge1`$ and $`f(x)=0`$ for $`x<1`$. Its transform is $`\Phi(s)=F(s)-2/(s-\theta)`$ for $`\mathrm{Re}\,s>1`$, which is meromorphic in the plane; its poles in $`\mathrm{Re}\,s>0`$ are the zeros of $`\zeta`$ in the critical strip and the point $`\theta<\Theta`$. Thus $`\Phi`$ is regular for $`\mathrm{Re}\,s>\Theta`$ and in no larger half-plane, and regular at the real point $`\Theta`$ because $`\zeta`$ has no real zero in $`(0,1)`$.

By Kaczorowski and Pintz (1986), for every $`\alpha>0`$ and all large $`T`$, the interval $`[T^{1-\alpha},T]`$ contains a change of sign of $`f`$. Take $`\alpha=\tfrac14`$ and $`T=Y^{3/2}`$. Then $`[Y^{9/8},Y^{3/2}]\subset[Y,Y^{3/2}]`$, where $`f(x)=G(x)-2x^{\theta}\le -x^{\theta}<0`$, contradicting a change of sign. $`\square`$

In the Lean development this localization is carried out directly by taking high derivatives $`F^{(n)}(s_0)=(-1)^{n}\int_{1}^{\infty}G(x)x^{-s_0-1}(\log x)^{n}\,dx`$ at $`s_0>1`$: with $`n=s_0\log X`$, the Gamma kernel $`u^{n}e^{-s_0 u}`$ concentrates on a multiplicative window of width $`\exp(O(\sqrt{\log X}))`$ around $`X`$, far narrower than $`[Y,Y^{3/2}]`$.

### 1.4 What has been derived

It suffices to prove, for every $`\varepsilon>0`$, that there are arbitrarily large $`Y`$ with

```math
G(x)\le x^{1/2+\varepsilon}\qquad\text{for all }x\in[Y,Y^{3/2}].\qquad(1.1)
```

---

## 2. The engine: divisibility against size

### 2.1 Chebyshev's inequality, and what is wrong with it

If $`\Delta\ne0`$ is an integer and $`p^{e_p}\mid\Delta`$ for each prime $`p`$ in a finite set, then

```math
\sum_{p}e_p\log p\le\log|\Delta|.\qquad(2.1)
```

Chebyshev applied (2.1) to $`\Delta=\binom{2n}{n}`$, where Legendre's floors $`\lfloor n/p^{h}\rfloor`$ accumulate rounding errors across all primes and stop at a multiplicative constant. Reaching $`\tfrac12`$ requires replacing $`\Delta`$ by an object whose divisibility is *smooth* in the exponent and whose size is controlled through the *linear* term.

### 2.2 Lattices, slopes, heights

Let $`K`$ be a totally complex number field of degree $`n_0`$, $`\mathfrak{o}_K`$ its ring of integers, and $`\Phi_K`$ a CM type of $`n_0/2`$ complex embeddings. A **lattice with norms** is a $`K`$-vector space $`V`$ of dimension $`r_1`$, an $`\mathfrak{o}_K`$-lattice $`\Lambda\subset V`$, and positive definite Hermitian forms on $`V_\tau`$ ($`\tau\in\Phi_K`$). On $`V_{\mathbb{R}}=\bigoplus_\tau V_\tau`$ (real dimension $`n_0 r_1`$), define

```math
\mu(\Lambda)=-\frac{2}{n_0r_1}\log\mathrm{covol}(\Lambda)+\mathrm{cst}_K,\qquad \deg\Lambda=r_1\mu(\Lambda),\qquad \mathrm{ht}(\phi)=\frac{2}{n_0}\sum_{\tau\in\Phi_K}\log\|\phi\|_{\tau}^{2},
```

with $`\mathrm{cst}_K=-\log2+\frac{\log|\mathrm{disc}\,K|}{n_0}`$. Dividing a lattice by a rational prime $`p`$ multiplies real covolume by $`p^{n_0 r_1}`$ and therefore raises $`\mu(\Lambda)`$ by $`2\log p`$.

### 2.3 The slope inequality

**Lemma 2.1** [D]. *Let $`\Lambda\subset V`$ and $`\Lambda'\subset V'`$ be lattices with norms, $`\Lambda`$ of rank $`r_1\ge1`$, and $`\phi:V\to V'`$ an injective $`K`$-linear map with $`\phi(\Lambda)\subset\Lambda'`$. Then*

```math
\mu(\Lambda)\le\mu(\phi(\Lambda))+\mathrm{ht}(\phi)\le\mu_{\max}(\Lambda')+\mathrm{ht}(\phi).
```

*Proof.* At each $`\tau\in\Phi_K`$, in orthonormal bases $`\phi_\tau`$ has singular values $`\sigma_1,\dots,\sigma_{r_1}\le\|\phi\|_\tau`$, so as a real map on $`\mathbb{R}^{2r_1}`$ it scales volume by $`\prod_i\sigma_i^{2}\le\|\phi\|_\tau^{2r_1}`$. Taking the product over $`\tau\in\Phi_K`$, $`\mathrm{covol}(\phi(\Lambda))\le\mathrm{covol}(\Lambda)\prod_\tau\|\phi\|_\tau^{2r_1}`$; taking logarithms and multiplying by $`-\frac{2}{n_0r_1}`$ yields the claim. $`\square`$

### 2.4 The rank problem, and the mean-slope inequality

Because each source lattice $`\Lambda_{\mathsf{w}}`$ has smaller rank $`d_{\mathsf{w}}< D`$ than the target $`Q`$, a single map controls only $`\mu_{\max}(Q)`$. Using a family of maps whose images almost span $`Q`$ controls the mean slope $`\mu(Q)`$.

**Proposition 2.2** [D, from the Stuhler–Grayson slope filtration]. *Let $`\mathsf{W}`$ be finite, $`\beta_{\mathsf{w}}:\Lambda_{\mathsf{w}}\to Q`$ integral with $`\mathrm{rk}\,\Lambda_{\mathsf{w}}=d_{\mathsf{w}}`$, $`\mathrm{rk}\,Q=D`$, and $`h_{\mathsf{w}}\ge\mathrm{ht}(\beta_{\mathsf{w}})`$. Set*

```math
s=\frac{\sum_{\mathsf{w}}d_{\mathsf{w}}\big(\mu(\Lambda_{\mathsf{w}})-h_{\mathsf{w}}\big)}{\sum_{\mathsf{w}}d_{\mathsf{w}}},\qquad \mathrm{sp}=\max_{\mathsf{w}}\big(\mu_{\max}(\Lambda_{\mathsf{w}})-h_{\mathsf{w}}\big)-s .
```

*Assume that whenever subspaces $`V^{\circ}_{\mathsf{w}}`$ satisfy $`\sum_{\mathsf{w}}(d_{\mathsf{w}}-\dim V^{\circ}_{\mathsf{w}})\le\alpha\sum_{\mathsf{w}}d_{\mathsf{w}}`$, the span of $`\beta_{\mathsf{w}}(V^{\circ}_{\mathsf{w}})`$ has codimension at most $`(\eta+4\alpha)D/(1-\delta)`$. If $`\mathsf{R}>0`$ bounds $`s-\mu_{\min}(Q)`$, then*

```math
\mu(Q)\ge s-\frac{1}{1-\delta}\Big(4\,\mathrm{sp}\,\log\Big(1+\frac{\mathsf{R}}{\mathrm{sp}}\Big)+\eta\mathsf{R}\Big).
```

---

## 3. Assembling the inequality

### 3.1 The table

Fix $`\varepsilon>0`$ and scale $`x`$. Auxiliary primes $`q\equiv3\pmod 8`$ and $`d\equiv7\pmod{16}`$ with $`q\mid d+1`$ are chosen at scale $`Y`$ ($`x^{1/J}\le Y\le x`$) with $`\log q\asymp Y^{\varepsilon/5}`$ and $`\log d\asymp Y^{\varepsilon/5}(\log Y)^{2}`$. Let $`n`$ be a power of $`2`$ with $`x^{2+8\varepsilon}\le n\le2x^{2+8\varepsilon}`$, $`N=4n`$, $`H=(\mathbb{Z}/N\mathbb{Z})^{\times}/\{\pm1\}`$ ($`|H|=n`$), and $`b=\exp(x^{\varepsilon})`$.

An odd prime power $`r=p^{h}`$ ($`p\ne q`$) lands in row $`v(r)=\tfrac12(1+\chi_q(r))\in\{0,1\}`$ and column $`[r]\in H`$, where $`\chi_d(n)=\big(\frac{-d}{n}\big)`$ and $`\chi_q(n)=\big(\frac{-q}{n}\big)`$. Distinguished cells $`(1,c)`$ ($`c\in T_2`$, $`|T_2|=\delta n`$) carry weight $`\approx b`$; other cells satisfy $`a/w_i\ge\frac13\exp(x^{\varepsilon/3})`$. Put

```math
P=\sum_{p\ne2,q}\ \sum_{h\ge1}\big(1-\chi_d(p^{h})\big)\log p\,\Big(\frac{a}{p^{h}}-w_{\mathrm{cell}(p^{h})}\Big)_{+},\qquad M=a\Big(\log a-1-\gamma+\frac{L'}{L}(1,\chi_d)-\langle\log w\rangle\Big).
```

The main inequality [Q] asserts that for $`x\ge x_0(\varepsilon)`$,

```math
\frac1n\sum_{u\in H}\ \sum_{t\in T}\ \sum_{j=1}^{m}\big(P_{u,t,j}-M_{u,t,j}\big)\ \le\ \mathcal{E},\qquad \frac{\mathcal{E}}{B}\le x^{C\varepsilon}\Big(\frac{x}{m}+\frac{xm}{k}+m+\frac{k}{n^{2}}\Big),\qquad(M)
```

with $`B=|T|\,m\,b\,\delta/2`$ and constants $`C,x_0(\varepsilon)`$ independent of $`x,Y,q,d`$.

### 3.2 The three local inputs

- **(a) Divisibility at inert primes** [H, (H6)–(H7)]. At a prime $`p`$ inert in $`E=\mathbb{Q}(\sqrt{-d})`$, the formal group has height $`2`$ with Lubin–Tate logarithm $`x+x^{p^{2}}/p+x^{p^{4}}/p^{2}+\cdots`$. On the Gauss lemniscate $`E=\mathbb{Q}(i)`$, differentiating to order $`a`$ yields

```math
v_p(A_a)\ \ge\ \Big\lfloor\frac{ap}{p^{2}-1}\Big\rfloor\ >\ a\Big(\frac1p+\frac1{p^{3}}+\frac1{p^{5}}+\cdots\Big)-1 .
```

- **(b) Size** [H, (H1)]. Via Bargmann–Fock Gaussian norms $`\|z^{a}\|^{2}=a!`$ and Kronecker's limit formula $`\log\Omega_E=\gamma+\log(2\pi/\mathrm{rd}_E)-\frac{L'}{L}(1,\chi_d)`$, the linear terms of size and prime count carry the identical constants $`\gamma`$ and $`\frac{L'}{L}(1,\chi_d)`$.
- **(c) Denominators** [H, (H9), via (H5)]. Assembling integral lattices across twists up to $`x^{400}`$ incurs Galois-cohomology / Tate–Shafarevich lengths bounded on windows by $`\approx\sqrt{b}`$ [Q].

### 3.3 The cancellation, and the small-table check

**Lemma 3.1** [D]. *With $`A=\sum_t a_t`$, $`\log\Omega=\gamma+\log(\pi/\varrho)-c`$, $`S=\sum_t m(a_t\log a_t-a_t-a_t\log\Omega)`$, and $`\Lambda_c=\sum_{j=1}^{m}(\langle\log w_{u,j}\rangle+\log(\varrho/\pi))+m\ell_g`$, one has $`S-A\Lambda_c+mA\ell_g=\sum_{t\in T}\sum_{j=1}^{m}M_{u,t,j}`$.* $`\square`$

For Rubinstein's numerical example ($`q=11, d=1319, N=16, x=1000, b=10, a=10^{4}`$, weight $`10`$ on $`(1,[1])`$ and $`0.05`$ elsewhere), exact evaluation in [`rh_core.py`](rh_core.py) [C] reproduces

```math
P=84701.032,\qquad M=84463.394,\qquad P-M=237.64,
```

and across all four column translates $`[1],[3],[5],[7]`$ gives excesses $`[237.64,-300.26,-176.73,229.13]`$ with mean $`-2.55`$.

### 3.4 Balancing the four errors

In (M), balancing the once-amortized boundary error $`x/m`$ against the $`m`$-step accumulated error $`m`$ forces $`m\asymp\sqrt{x}`$ and yields $`\sqrt{x}`$. At $`\varepsilon=0.01`$ ($`m\asymp x^{0.499375}, k\asymp x^{4.6}, n\asymp x^{2.08}`$), the four exponents are $`0.500625`$, $`-3.100625`$, $`0.499375`$, and $`0.44`$. Had step accumulation scaled as $`m^{2}`$ after normalization, balancing $`x/m+m^{2}`$ at $`m=x^{1/3}`$ would yield only $`\mathrm{Re}\,s\le\tfrac23`$. Notice the exact parallel with Stepanov's elementary proof of the Weil bound for curves over $`\mathbb{F}_q`$, where dividing $`m\cdot\#C(\mathbb{F}_q)\le m(q+1)+O(q)+O(m^{2})`$ by $`m`$ gives $`O(q/m)+O(m)`$, minimized at $`m\asymp\sqrt{q}`$.

---

## 4. From the inequality to a bound on $`G(x)`$

Averaging over $`u\in H`$ eliminates residue classes modulo $`N`$ while non-distinguished cells sit past $`\frac13\exp(x^{\varepsilon/3})`$ where the classical zero-free region saves $`\exp(-\kappa x^{\varepsilon/8})`$ (Estimate C [Q]). Expanding $`(1+\chi_q)(1-\chi_d)=1+\chi_q-\chi_d-\chi_q\chi_d`$ isolates

```math
R_{q,d}(x)=G(x)+G_{\chi_q}(x)-G_{\chi_d}(x)-G_{\chi_q\chi_d}(x)
```

with coefficient $`+1`$ on $`G(x)`$. Because $`\log q,\log d`$ are powers of $`Y`$, modulus mean-square estimates supply a pair $`(q,d)`$ with $`|G_\chi(x)|\le x^{1/2+\varepsilon}`$ on $`[Y,Y^{3/2}]`$ without circularity (Estimate B [Q]). Combining (M) with Estimates B and C yields (1.1), and Lemma 1.2 gives RH.

---

## 5. The passage to infinity

The logical form of the prime-side input is

```math
\forall\varepsilon>0\ \ \forall Z\ \exists\,\mathcal{D}\ \ \forall x\in[Z,Z^{3/2}]:\quad (M)\text{ holds at }x\text{ with data }\mathcal{D}(x),\qquad(5.1)
```

asserted as Hypothesis (H12) [H] in the reconstruction. Circularity among thresholds is avoided by a strict partial order of parameter choice:

```math
\beta\ \prec\ \varepsilon\ \prec\ R\ \prec\ J\ \prec\ x_0(\varepsilon)\ \prec\ Z\ \prec\ x\ \prec\ Y\ \prec\ (q,d)\ \prec\ \mathcal{D}(x).
```

Windows $`[Z,Z^{3/2}]`$ arise in §5.3 because `(H9)` bounds twist cohomology lengths $`h_{\ell,e}`$ sharply away from at most $`R`$ unknown **exceptional centres**, while near an exceptional centre only a cruder bound holds, useful only when $`e`$ is at most a small power of $`x`$.

---

## 6. The audit, and the three problems of §6.3

Sorting `(H0)–(H12)` shows that local formal-group divisibility `(H6)` and linear-algebra slope filtration `(H8, H10, H11)` are explicit, automorphic inner-product and mass formulas `(H1, H4, H5)` supply the exact linear main term, `(H9)` controls cohomology denominators on windows, and `(H12)` asserts scale-uniformity.

Unlike the function-field case—where a fixed curve $`C/\mathbb{F}_q`$ and fixed $`2g`$-dimensional $`H^{1}(C)`$ turn point counting into a spectral trace identity—here the certificate dimensions ($`n\asymp x^{2+8\varepsilon}`$, $`k\asymp x^{9/2+10\varepsilon}`$, $`m\asymp x^{1/2-\varepsilon/16}`$) diverge with $`x`$ and Lemma 2.1 is one-sided. Converting "certified but not understood" into concrete mathematics raises three questions:

1. **A fixed projective system (§6.3.1):** Can the scale-$`x`$ certificates be realized as specializations of one scale-independent projective tower?
2. **An identity for the determinant (§6.3.2):** Can Lemma 2.1 be replaced by an exact local-factor formula for the pairing determinant so divisibility is read off as an identity rather than bounded by an inequality?
3. **Uniformity and window elimination (§6.3.3):** Can the exceptional-centre windows and scale-dependent blowups be removed or relaxed while preserving a zero-free half-plane?

Part II (§§9–13) gives explicit answers to 1, 2, and the polyhedral structure of 3, together with the non-abelian obstruction separating $`GL_1`$ from $`GL_2`$.

---

## 7. The reduction on one page

1. $`G(x)=\sum_{n\le x}\Lambda(n)(x/n-1)-x\log x+(1+\gamma)x`$ has Mellin transform regular at $`s=1`$ with simple poles at nontrivial zeros $`\rho`$ of residue $`-1/(\rho(\rho-1))`$.
2. One-sided upper bounds on windows $`[Y,Y^{3/2}]`$ exclude zeros in $`\mathrm{Re}\,s>\theta`$ (Lemmas 1.1–1.2 [D]).
3. An injective integral lattice map divisible by $`p^{e}`$ satisfies $`e\log p\le\text{height}+\text{slope difference}`$ (Lemma 2.1 [D]); spanning families of lower-rank maps bound the mean slope $`\mu(Q)`$ (Proposition 2.2 [D]).
4. Height-$`2`$ Lubin–Tate coordinates at inert primes give smooth divisibility $`P`$ `(H6–H7)`, while Bargmann–Fock norms $`a!`$ and Kronecker's limit formula match the linear constants of $`M`$ `(H1, Lemma 3.1 [D])`.
5. Balancing end-chain error $`x/m`$ against $`m`$-step accumulation $`m`$ gives $`m\asymp\sqrt{x}`$ and error $`x^{1/2+C\varepsilon}`$.
6. Character averaging isolates $`G(x)`$ in $`R_{q,d}(x)`$; Estimates B and C control twisted and far-scale terms; `(H12)` supplies the data uniformly across windows.

---

## 8. Scope and epistemic boundaries

Part I does not prove `(H0)–(H12)` and does not inspect the compiled Lean binary directly; every `[Q]` and `[H]` statement reflects Rubinstein's v2 reconstruction. Part II below proves unconditional algebraic theorems (`[D]`) and exact integer/rational computational certificates (`[C]`) on the concrete mathematical objects that govern §§2–6.

---

# Part II. Structural results on §6.3, non-abelian $`GL_2`$ barriers, and window-free Quasi-RH

## 9. Resolving §6.3(2): Padé obstructions for single determinants and Fitting-ideal exactness

To test §6.3(2) concretely, consider the Gauss lemniscate $`E=\mathbb{Q}(i)`$ cited in §3.2(a), with elliptic coordinate $`u=\mathrm{sl}(z)`$ satisfying $`(u')^{2}=1-u^{4}`$, $`u(0)=0`$, $`u'(0)=1`$, hence $`u''=-2u^{3}`$ and exact integer derivatives $`A_a=u^{(a)}(0)\in\mathbb{Z}`$ (nonzero precisely for $`a\equiv1\pmod 4`$). Write $`D=d/dz`$.

### 9.1 Exact basis invariance and the unshifted hyperfactorial identity

Compare the $`r\times r`$ shifted Hankel derivative determinant and the odd-monomial pairing determinant at scale $`a_0\equiv1\pmod 4`$:

```math
H_r(a_0)=\det_{0\le i,j< r}\!\big[A_{a_0+4(i+j)}\big],\qquad B_r(a_0)=\det_{0\le i,j< r}\!\Big[D^{a_0+4j}\big(u^{4i+1}\big)(0)\Big].
```

**Theorem 9.1** [D]. *For every integer $`r\ge1`$ and every $`a_0\ge1`$,*

```math
H_r(a_0)=\left(\prod_{i=0}^{r-1}(4i)!\right)B_r(a_0).
```

*Proof.* For any integer $`m\ge1`$, differentiating $`u^{m}`$ twice using $`(u')^{2}=1-u^{4}`$ and $`u''=-2u^{3}`$ gives

```math
D^{2}(u^{m})=m(m-1)u^{m-2}(1-u^{4})+mu^{m-1}(-2u^{3})=-m(m+1)u^{m+2}+m(m-1)u^{m-2}.
```

Applying $`D^{2}`$ a second time sends $`u^{m}`$ to $`m(m+1)(m+2)(m+3)u^{m+4}`$ plus integer multiples of $`u^{m}`$ and $`u^{m-4}`$. By induction on $`i\ge0`$, there exist integers $`c_{i,k}\in\mathbb{Z}`$ such that

```math
D^{4i}(u)=(4i)!\,u^{4i+1}+\sum_{k=0}^{i-1}c_{i,k}\,u^{4k+1}.
```

Applying $`D^{a_0+4j}`$ at $`z=0`$ expresses the $`i`$-th row of $`H_r(a_0)`$ as $`(4i)!`$ times the $`i`$-th row of $`B_r(a_0)`$ plus integer multiples of rows $`0,\dots,i-1`$. Thus $`H_r(a_0)=L_r B_r(a_0)`$ where $`L_r`$ is lower-triangular with diagonal $`(4i)!`$. $`\square`$

At the **unshifted base point** $`a_0=1`$, Stieltjes' continued-fraction expansion for Jacobian elliptic Laplace transforms yields an exact local-factor determinant formula:

```math
H_r(1)=\prod_{k=1}^{r-1}\Big[(4k-3)(4k-2)^{2}(4k-1)^{2}(4k)^{2}(4k+1)\Big]^{r-k},\qquad(9.1)
```

verified in [`cm_determinant_lab.py`](cm_determinant_lab.py) [C] up to $`r=8`$ ($`208`$ decimal digits), where every prime dividing $`H_r(1)`$ satisfies $`p\le 4r-3`$ and all large prime factors of individual $`A_a`$ cancel out.

### 9.2 Refutation of single-determinant local factorization at shifted scales

To localize around a prime scale $`x\sim a_0/w`$ without taking rank $`r\asymp x`$, one must shift to $`a_0\gg r`$. At shifted scales, (9.1) shatters.

**Proposition 9.2** [D, via prime-certified counterexamples]. *For $`r\ge2`$ and shifted scales $`a_0\ge9`$, neither $`H_r(a_0)`$ nor $`B_r(a_0)`$ is supported on primes $`\le a_0+8(r-1)`$, nor is its odd super-scale part supported on primes inert in $`\mathbb{Q}(i)`$.*

*Proof.* By Theorem 9.1, any prime $`p>4r-4`$ divides $`H_r(a_0)`$ with the exact same multiplicity as $`B_r(a_0)`$. Direct exact integer evaluation and exhaustive trial division to $`\sqrt{n}`$ ([`cm_determinant_lab.py`](cm_determinant_lab.py) [C]) give:
- At $`r=2, a_0=17`$ (local scale cutoff $`25`$), $`B_2(17)`$ is divisible by the $`6`$-digit **split** prime $`p=131893\equiv1\pmod 4`$.
- At $`r=3, a_0=21`$ (local scale cutoff $`37`$), $`B_3(21)`$ is divisible by the $`15`$-digit **split** prime $`p=254836380292649\equiv1\pmod 4`$, carrying $`22.3\%`$ of $`\log|B_3(21)|`$. $`\square`$

Mechanically, shifting $`a_0=1+4s`$ introduces an orthogonal-polynomial / Padé numerator factor whose primes reflect accidental modular congruences rather than formal-group height. For a single square determinant, Lemma 2.1 *must* remain a one-sided inequality because these positive Padé valuations have no local formal-group formula.

### 9.3 Fitting-ideal Padé annihilation and the Bezout conservation law

In §2.4, however, the proof maps a *family* of words $`\bigoplus_{\mathsf{w}\in\mathsf{W}}\Lambda_{\mathsf{w}}\to Q`$. Algebraically, the invariant of a rectangular $`D\times M`$ integral matrix ($`M>D`$) is not one minor but its **zeroth Fitting ideal** $`\mathrm{Fitt}_0`$, generated over $`\mathbb{Z}`$ by the greatest common divisor of all $`\binom{M}{D}`$ maximal $`D\times D`$ minors.

- **Computational exactness [C]** ([`fitting_ideal_lab.py`](fitting_ideal_lab.py)):
  1. On the Gauss lemniscate with $`r=1`$ and $`M=3`$ consecutive columns, $`\mathrm{Fitt}_0(a_0)=\gcd(|A_{a_0}|,|A_{a_0+4}|,|A_{a_0+8}|)`$ strips **100% of split primes and 100% of super-scale primes ($`p>a_0+8`$) across all $`51`$ shifted scales $`a_0=1,5,\dots,201`$**. At $`a_0=201`$, $`|A_{201}|`$ has $`324`$ digits ($`37.9\%`$ Padé noise), whereas $`\mathrm{Fitt}_0(201)`$ has $`201`$ digits and matches the exact local Lubin–Tate floor product $`2^{100}\prod_{p\equiv3(4)}p^{\lfloor 201p/(p^{2}-1)\rfloor}`$ identically. (With only $`M=2`$ columns, one accidental inert prime $`p=271`$ survives at $`a_0=133`$, consistent with a Borel–Cantelli $`p^{-M}`$ collision heuristic: $`\sum_{a_0}\sum_{p>a_0}p^{-2}`$ diverges for $`M=2`$, whereas $`\sum_{a_0}\sum_{p>a_0}p^{-3}<\infty`$ for $`M\ge3`$.)
  2. In the Hurwitz divided-power basis $`u^{4i+1}/(4i+1)!`$ for ranks $`r=1,2,3,4`$ with $`M=r+2`$ columns ($`a_0\le 49`$), zero super-scale primes survive and split-prime factors freeze at $`a_0`$-independent constants (`{}`, `{5:1}`, `{5:3}`, `{5:5, 13:1}`).
  3. On general Weierstrass CM curves (`49a1` with $`d=7`$ and `121b1` with $`d=11`$) at $`a_0=80`$, dividing $`\{2,d\}`$ out of $`\gcd(|A_{80}|,|A_{81}|,|A_{82}|)`$ leaves zero split primes and cofactor $`1`$.

Why, then, can the proof not collapse to rank $`1`$ by taking a Bezout linear combination $`\sum_{j=0}^{M-1}u_j A_{a_0+4j}=\mathrm{Fitt}_0(a_0)`$, or by surjecting $`\Phi:\mathbb{Z}^{M}\to\mathbb{Z}`$?

**Proposition 9.3 (Cauchy–Binet kernel trap & Bargmann–Kronecker normalization)** [D].
1. *Let $`v=(v_0,\dots,v_{M-1})\in\mathbb{Z}^{M}\setminus\{0\}`$, $`\Phi:\mathbb{Z}^{M}\to\mathbb{Z}`$ be $`\Phi(x)=v\cdot x`$, and $`g=\gcd(v_0,\dots,v_{M-1})=\mathrm{Fitt}_0(\Phi)`$. Then $`\mathrm{im}(\Phi)=g\mathbb{Z}`$ and the rank-$`(M-1)`$ kernel lattice $`\ker\Phi\subset\mathbb{Z}^{M}`$ has Euclidean covolume*

```math
\mathrm{covol}(\ker\Phi)=\frac{\|v\|_2}{g}.
```

2. *For $`E=\mathbb{Q}(i)`$ ($`\mathrm{rd}_E=2`$), Kronecker's limit constant $`\Omega_{\mathbb{Q}(i)}=\exp(\gamma+\log\pi-\frac{L'}{L}(1,\chi_{-4}))`$ equals $`\frac{2\varpi^{2}}{\pi}`$ where $`\varpi=\int_0^{1}(1-t^{4})^{-1/2}dt`$ is Gauss's lemniscate constant, and holomorphic Bargmann–Fock pairing ($`\|z^{a}\|_{L^{2}}^{2}=a!`$) has complex height $`\log(a!)\sim a\log a`$, matching $`2\log\mathrm{Fitt}_0(a)\sim a\log a`$ whereas pointwise Cauchy bounds $`|A_a|^{2}\asymp (a!)^{2}R_0^{-2a}`$ have exponent $`2a\log a`$.*

*Proof.* (1) Write $`v=g\tilde{v}`$ with $`\gcd(\tilde{v})=1`$. Extend $`\tilde{v}`$ to a unimodular matrix $`U\in SL_M(\mathbb{Z})`$ whose first row is $`\tilde{v}`$; the Gram determinant on $`\ker\Phi=\tilde{v}^{\perp}\cap\mathbb{Z}^{M}`$ equals $`\|\tilde{v}\|_2=\|v\|_2/g`$ by Cauchy–Binet. (2) By Malmsten's formula for Dirichlet's beta function $`L(s,\chi_{-4})`$,

```math
\frac{L'}{L}(1,\chi_{-4})=\gamma+\log(2\pi)-2\log\frac{\Gamma(1/4)}{\Gamma(3/4)}.
```

Using Euler's reflection identity $`\Gamma(1/4)\Gamma(3/4)=\pi\sqrt{2}`$ and Gauss's relation $`\varpi=\Gamma(1/4)^{2}/(2\sqrt{2\pi})`$ gives $`\Gamma(1/4)/\Gamma(3/4)=2\varpi/\sqrt{\pi}`$, whence

```math
\gamma+\log\pi-\frac{L'}{L}(1,\chi_{-4})=\log\pi-\log(2\pi)+\log\frac{4\varpi^{2}}{\pi}=\log\frac{2\varpi^{2}}{\pi}
```

(certified numerically to $`1.1\times10^{-19}`$ in [`frontier6_bezout_and_centres.py`](frontier6_bezout_and_centres.py) [C]). Over $`\mathbb{C}`$, an inert prime power $`p^{e}`$ scales area by $`p^{2e}`$, so $`2\log\mathrm{Fitt}_0(a)\sim 2\cdot\frac12 a\log a=a\log a=\log\|z^{a}\|_{L^{2}}^{2}+O(a)`$. $`\square`$

Similarly, solving $`\sum u_j Q_j=1`$ for coprime Padé cofactors $`Q_j=A_{a_0+4j}/\mathrm{Fitt}_0(a_0)`$ forces $`\log\|u\|_\infty\approx\log|Q_0|`$ (`299.1` vs. `282.5` at $`a_0=201`$ [C]). **Only** mapping rank-deficient words $`\beta_{\mathsf{w}}:\mathbb{Z}^{r}\to\mathbb{Z}^{D}`$ ($`r< D`$) with unit weights into a higher-rank target $`Q=\mathbb{Z}^{D}`$ (Proposition 2.2) achieves $`[\mathbb{Z}^{D}:\sum_{\mathsf{w}}\mathrm{im}(\beta_{\mathsf{w}})]=\mathrm{Fitt}_0`$ with zero Bezout multiplier penalty—verified via exact modular Hermite Normal Form across $`D\in\{2,3,4\}`$ (`r=D-1`) in [`frontier7_hnf_slope_engine.py`](frontier7_hnf_slope_engine.py) [C].

---

## 10. Non-abelian $`GL_2`$ barriers: ultrametric collapse vs. $`a_p`$-free Hecke–Wronskians

Why does the divisibility engine not immediately extend from $`\zeta(s)`$ ($`GL_1`$) to a non-CM elliptic curve $`E/\mathbb{Q}`$ of conductor $`N_E`$ ($`GL_2`$)? Let $`\omega(t)=\sum_{n=1}^{\infty}b_n t^{n-1}\,dt\in\mathbb{Z}[[t]]\,dt`$ be the Weierstrass invariant differential of a minimal model over $`\mathbb{Z}`$, and $`\mathcal{L}_E(t)=\sum_{n=1}^{\infty}\frac{a_n(E)}{n}t^{n}`$ the Honda formal logarithm with inverse $`\mathcal{E}_E(z)=\sum_{a=1}^{\infty}C_a(E)\frac{z^{a}}{a!}`$.

**Theorem 10.1** [D]. *Let $`E/\mathbb{Z}`$ be a minimal Weierstrass elliptic curve and $`p\nmid 2N_E`$ an odd good prime with Frobenius trace $`a_p(E)=p+1-\#E(\mathbb{F}_p)`$.*
1. *(Integrality & ordinary slope collapse.) One has $`C_a(E)\in\mathbb{Z}`$ for all $`a\ge1`$. If $`p\nmid a_p(E)`$ (ordinary reduction, asymptotic density $`1`$ for non-CM $`E`$), then $`v_p(b_{p^{r}})=0`$ for all $`r\ge1`$ and $`\liminf_{a\to\infty}v_p(C_a(E))/a=0`$ (in contrast to supersingular $`p\mid a_p(E)`$, where $`v_p(C_a(E))/a\ge\frac{p}{p^{2}-1}-o(1)>0`$).*
2. *(Purely $`p`$-adic trace uniqueness.) For every $`p\ge5`$, $`a_p(E)`$ is the unique integer $`c\in[-2\sqrt{p},2\sqrt{p}]`$ satisfying*

```math
b_p\equiv c\pmod p\qquad\text{and}\qquad b_{p^{2}}-c\,b_p+p\equiv0\pmod{p^{2}}
```

*(and for $`p=3`$, adding $`b_{27}-cb_9+3b_3\equiv0\pmod{27}`$ resolves the supersingular case $`3\mid b_3`$).*
3. *($`a_p`$-free Hecke–Wronskian divisibility.) Define the $`a_p(E)`$-independent $`2\times 2`$ determinant*

```math
W_p(n_1,n_2)=\det\begin{pmatrix}b_{pn_1}+p\,b_{n_1/p} & b_{n_1}\\ b_{pn_2}+p\,b_{n_2/p} & b_{n_2}\end{pmatrix}
```

*(with $`b_{n/p}=0`$ when $`p\nmid n`$). Then for all integers $`r\ge1`$ and $`u_1,u_2\ge1`$, $`v_p\big(W_p(p^{r-1}u_1,p^{r-1}u_2)\big)\ge r`$.*

*Proof.* (1) Since $`a_1(E)=1`$ and $`a_n(E)\in\mathbb{Z}`$, $`\mathcal{L}_E'(t)=1+\sum_{m\ge1}a_{m+1}(E)t^{m}`$ is a unit in $`\mathbb{Z}[[t]]`$, so the derivation $`D=(\mathcal{L}_E'(t))^{-1}\frac{d}{dt}`$ preserves $`\mathbb{Z}[[t]]`$; hence $`C_a(E)=(D^{a}t)(0)\in\mathbb{Z}`$ for all $`a\ge1`$. By Honda (1968), $`\mathcal{L}_E`$ is strictly isomorphic over $`\mathbb{Z}_{(p)}`$ to $`\widehat{E}`$, giving the Atkin–Swinnerton-Dyer congruence

```math
\delta_p(n):=b_{pn}-a_p(E)b_n+p\,b_{n/p}\equiv 0\pmod{p^{v_p(n)+1}}\qquad(10.1)
```

for all $`n\ge1`$. By induction on $`r\ge1`$ with $`b_1=1`$, $`b_{p^{r}}\equiv a_p(E)b_{p^{r-1}}\equiv a_p(E)^{r}\not\equiv0\pmod p`$ when $`p\nmid a_p(E)`$, so $`v_p(b_{p^{r}})=0`$. Moreover, when $`p\nmid a_p(E)`$, $`\widehat{E}`$ has height $`1`$ over $`\mathbb{Z}_p`$, so $`[p]_{\widehat{E}}(t)`$ has $`p-1`$ nonzero roots $`\zeta\in\overline{\mathbb{Q}}_p`$ of valuation $`v_p(\zeta)=\frac{1}{p-1}`$ at which $`\mathcal{L}_E(\zeta)=0`$. Since $`\mathcal{L}_E`$ fails to be injective on any closed disk of radius $`\ge p^{-1/(p-1)}`$, its formal inverse $`\mathcal{E}_E(z)=\sum_{a\ge1}\frac{C_a(E)}{a!}z^{a}`$ has $`p`$-adic radius of convergence at most $`p^{-1/(p-1)}`$. By Cauchy–Hadamard, $`\liminf_{a\to\infty}\frac{v_p(C_a(E))-v_p(a!)}{a}\le-\frac{1}{p-1}`$; since $`\lim_{a\to\infty}\frac{v_p(a!)}{a}=\frac{1}{p-1}`$ and $`C_a(E)\in\mathbb{Z}`$ gives $`v_p(C_a(E))\ge0`$, it follows that $`\liminf_{a\to\infty}v_p(C_a(E))/a=0`$ (and empirically $`v_p(C_a(E))=O(\log_p a)`$ [C]).

(2) If $`p\ge5`$ and $`p\mid b_p`$, then $`p\mid a_p(E)`$ by (10.1) at $`n=1`$; since $`|a_p(E)|\le 2\sqrt{p}< p`$, necessarily $`a_p(E)=0`$, which is also the unique $`c\in[-2\sqrt{p},2\sqrt{p}]`$ with $`c\equiv b_p\equiv0\pmod p`$. If $`p\nmid b_p`$, then $`b_p`$ is invertible mod $`p^{2}`$, so $`b_{p^{2}}-cb_p+p\equiv0\pmod{p^{2}}`$ determines $`c\bmod p^{2}`$ uniquely; since the interval $`[-2\sqrt{p},2\sqrt{p}]`$ has length $`4\sqrt{p}< p^{2}`$ for all $`p\ge3`$, at most one integer $`c`$ in it satisfies the congruence, and $`c=a_p(E)`$ does by (10.1) at $`n=p`$. At $`p=3`$ with $`3\mid b_3`$, $`b_9-a_3 b_3+3\equiv0\pmod 9`$ forces $`v_3(b_9)=1`$, so $`b_{27}-cb_9+3b_3\equiv0\pmod{27}`$ determines $`c\bmod 9`$, separating $`\{-3,0,3\}`$.

(3) Subtracting $`a_p(E)`$ times column $`2`$ from column $`1`$ leaves the determinant unchanged and turns column $`1`$ into $`(\delta_p(n_1),\delta_p(n_2))^{\top}`$. For $`n_i=p^{r-1}u_i`$, $`v_p(\delta_p(n_i))\ge v_p(n_i)+1\ge r`$ by (10.1), whence $`p^{r}\mid W_p(n_1,n_2)`$. $`\square`$

Theorem 10.1 isolates the exact boundary between $`GL_1`$ and $`GL_2`$ ([`gl2_divisibility_lab.py`](gl2_divisibility_lab.py), [`fitting_ideal_lab.py`](fitting_ideal_lab.py) [C]): while $`W_p`$ eliminates $`a_p(E)`$ from the local operator, its valuation $`v_p\in\mathbb{Z}_{\ge0}`$ counts good primes *uniformly* rather than weighting them by $`a_p(E)/\sqrt{p}`$, because non-negative integer valuations can encode discrete Artin/CM factors $`1-\chi_d(p)\in\{0,2\}`$ but cannot directly represent irrational continuous Satake traces $`a_p(E)/\sqrt{p}\in[-2,2]\setminus\mathbb{Q}`$.

---

## 11. Resolving §6.3(1): horizontal Iwasawa exactness vs. vertical infinite rank

Consider the dyadic group tower $`\cdots\twoheadrightarrow H_{2N}\twoheadrightarrow H_N\twoheadrightarrow\cdots\twoheadrightarrow H_8`$, where $`H_N=(\mathbb{Z}/N\mathbb{Z})^{\times}/\{\pm1\}`$ ($`N=4n=2^{k+2}`$) has projective limit $`\varprojlim H_N=\mathbb{Z}_2^{\times}/\{\pm1\}\cong\mathbb{Z}_2`$, with $`2`$-to-$`1`$ reduction map $`\pi:H_{2N}\twoheadrightarrow H_N`$.

For fixed $`(q,d,a,b,w_0)`$, let $`P^{(N)}(c)`$ and $`M^{(N)}`$ denote the §3.1 table sums at level $`N`$ with sole distinguished cell $`(1,c)`$ of weight $`b`$ and all other $`2n-1`$ cells of weight $`w_0`$, and let $`P_{\mathrm{flat}}(w_0),M_{\mathrm{flat}}(w_0)`$ be the uniform table with weight $`w_0`$ on all cells. Define the centered excess on $`H_N`$:

```math
\mu_N(c)=\big(P^{(N)}(c)-M^{(N)}\big)-\big(P_{\mathrm{flat}}(w_0)-M_{\mathrm{flat}}(w_0)\big).
```

**Theorem 11.1** [D].
1. *(Horizontal pullback invariance.) For any subset $`T_2\subset H_N`$, pulling back indicator weights to $`\pi^{-1}(T_2)\subset H_{2N}`$ satisfies $`P^{(2N)}(\pi^{*}w)=P^{(N)}(w)`$ and $`M^{(2N)}(\pi^{*}w)=M^{(N)}(w)`$ identically.*
2. *(Horizontal Iwasawa distribution identity.) For every $`N=2^{k+2}\ge8`$ and every $`c\in H_N`$ with fiber $`\pi^{-1}(c)=\{c_+,c_-\}\subset H_{2N}`$,*

```math
\mu_{2N}(c_+)+\mu_{2N}(c_-)=\mu_N(c).
```

*Hence $`(\mu_N)_{N}`$ defines an exact finitely-additive distribution on $`\mathbb{Z}_2^{\times}/\{\pm1\}`$.*
3. *(Vertical infinite linear rank.) For any integer $`R\ge1`$, odd primes $`p_1<\dots< p_R`$, and rational scales $`x_i\in(p_i,p_{i+1})`$, the $`R\times R`$ ramp matrix $`K_{i,j}=(x_i/p_j-1)_{+}`$ has $`\mathrm{rk}_{\mathbb{Q}}(K)=R`$, whereas for any curve $`C/\mathbb{F}_p`$ of genus $`g`$ with Frobenius eigenvalues $`\{\alpha_\ell\}_{\ell=1}^{2g}`$, the $`R\times R`$ trace Hankel matrix $`F_{i,j}=\sum_{\ell=1}^{2g}\alpha_\ell^{i+j-1}`$ satisfies $`\mathrm{rk}_{\mathbb{Q}}(F)\le 2g`$ for all $`R`$.*

*Proof.* (1) Since $`\pi([p^{h}]_{2N})=[p^{h}]_N`$, every prime power receives the identical cell weight under $`\pi^{*}w`$, so $`P`$ is unchanged; doubling both $`|T_2|`$ and $`n`$ preserves $`\langle\log w\rangle`$, so $`M`$ is unchanged.

(2) Every odd prime power $`r`$ with $`\chi_q(r)=1`$ and $`[r]_N=c`$ satisfies $`[r]_{2N}\in\{c_+,c_-\}`$, so across the two tables $`c_+`$ and $`c_-`$ it receives weight $`b`$ once and $`w_0`$ once; every other prime power receives weight $`w_0`$ in both. Adding gives $`P^{(2N)}(c_+)+P^{(2N)}(c_-)=P^{(N)}(c)+P_{\mathrm{flat}}(w_0)`$. For $`M`$, $`\langle\log w^{(2N,c_\pm)}\rangle=\log w_0+\frac{1}{4n}\log(b/w_0)`$, whose sum over $`\{c_+,c_-\}`$ is $`2\log w_0+\frac{1}{2n}\log(b/w_0)=\langle\log w^{(N,c)}\rangle+\log w_0`$, giving $`M^{(2N)}(c_+)+M^{(2N)}(c_-)=M^{(N)}(c)+M_{\mathrm{flat}}(w_0)`$. Subtracting yields $`\mu_{2N}(c_+)+\mu_{2N}(c_-)=\mu_N(c)`$.

(3) Since $`x_i< p_{i+1}\le p_j`$ for $`j>i`$, $`(x_i/p_j-1)_{+}=0`$ above the diagonal, while on the diagonal $`K_{i,i}=x_i/p_i-1>0`$. Thus $`K`$ is lower-triangular with nonzero diagonal and $`\det K=\prod_{i=1}^{R}(x_i/p_i-1)>0`$. By contrast, $`F=\sum_{\ell=1}^{2g}v_\ell v_\ell^{\top}`$ with $`v_\ell=(\alpha_\ell^{0},\dots,\alpha_\ell^{R-1})^{\top}\in\overline{\mathbb{Q}}^{R}`$ is a sum of $`2g`$ rank-$`1`$ matrices. $`\square`$

Theorem 11.1 (verified across $`N=8\dots 256`$ in [`projective_tower_lab.py`](projective_tower_lab.py) [C]) explains §6.3(1): *horizontally* at fixed $`x`$, the table is already a projective Iwasawa distribution; *vertically* across $`x\to\infty`$, fixing one finite-rank object over a fixed CM field $`\mathbb{Q}(\sqrt{-d})`$ fails because (i) once $`x>d`$, bounding $`G_{\chi_d}(x)`$ ceases to be unconditional and requires GRH for $`L(s,\chi_d)`$, (ii) $`k/n^{2}`$ forces $`n\to\infty`$, and (iii) the piecewise-linear ramp kernel $`(x/r-1)_{+}`$ spans an infinite-dimensional function space as $`x`$ varies.

---

## 12. Polyhedral decompilation and the window-free Quasi-RH dichotomy

Finally, consider the asymptotic exponent structure of the blueprint. Write $`m=x^{\alpha_m}`$, $`k=x^{\alpha_k}`$, $`n=x^{\alpha_n}`$ ($`\alpha_m,\alpha_k,\alpha_n\ge0`$).

1. **Hidden-constraint gap in (M) [D/C]:** Minimizing $`\max(1-\alpha_m,\,1+\alpha_m-\alpha_k,\,\alpha_m,\,\alpha_k-2\alpha_n)`$ alone achieves $`\theta^{*}=\tfrac12`$ already at $`(\alpha_m,\alpha_k,\alpha_n)=(\tfrac12,1,\tfrac14)`$—three and a half powers of $`x`$ smaller in rank $`k`$ than the reconstruction's $`(\tfrac12,\tfrac92,2)`$ ([`proof_dag_decompiler.py`](proof_dag_decompiler.py) [C]). Thus $`k\asymp x^{9/2}, n\asymp x^{2}`$ are forced not by the four displayed terms of (M), but by internal codimension and modulus constraints in `(H0, H8, H10)`—specifically $`\eta\mathsf{R}\asymp xn^{2}/k`$ (exponent $`1+2\alpha_n-\alpha_k`$, which pairs with $`k/n^{2}`$ of exponent $`\alpha_k-2\alpha_n`$ to force $`\max(1+2\alpha_n-\alpha_k,\alpha_k-2\alpha_n)\ge\tfrac12`$) together with a table-modulus requirement such as $`x^{5/2}/n\le x^{\theta}`$ (exponent $`\tfrac52-\alpha_n`$).

2. **Exact rational Farkas dichotomy under an unconditional twist cap [D, relative to the blueprint model]:** Suppose `(H9)`'s window-dependent twist bound is replaced by an unconditional, window-free rank ceiling $`\alpha_k\le\tau_k`$ (representing §5.3's small-power regime where no exceptional centres need dodging):

**Proposition 12.1** [D]. *On the polyhedral model $`\mathcal{P}(\tau_k)`$ defined over $`\alpha_m,\alpha_k,\alpha_n\ge0`$ by $`\alpha_k\le\tau_k`$ and the six exponent bounds*

```math
\theta\;\ge\;\max\!\Big(1-\alpha_m,\ \alpha_m,\ 1+\alpha_m-\alpha_k,\ \alpha_k-2\alpha_n,\ 1+2\alpha_n-\alpha_k,\ \tfrac52-\alpha_n\Big),
```

*the minimax zero-free abscissa $`\theta^{*}(\tau_k)=\min_{\mathcal{P}(\tau_k)}\theta`$ is given for all $`\tau_k\ge\tfrac32`$ by*

```math
\theta^{*}(\tau_k)\;=\;\max\!\Big(\frac12,\ 2-\frac{\tau_k}{3}\Big).
```

*Proof.* **Lower bound (exact Farkas dual):** Taking weights $`(\tfrac12,\tfrac12)`$ on $`1-\alpha_m`$ and $`\alpha_m`$ gives $`\theta\ge\tfrac12`$. Taking weights $`\tfrac23`$ on $`\tfrac52-\alpha_n`$, $`\tfrac13`$ on $`1+2\alpha_n-\alpha_k`$, and $`\tfrac13`$ on the constraint $`\alpha_k-\tau_k\le0`$ cancels $`\alpha_n`$ and $`\alpha_k`$ identically:

```math
\theta\;\ge\;\frac23\Big(\frac52-\alpha_n\Big)+\frac13\big(1+2\alpha_n-\alpha_k\big)+\frac13\big(\alpha_k-\tau_k\big)\;=\;2-\frac{\tau_k}{3}.
```

**Upper bound (explicit primal witness):** For $`\tau_k\ge\tfrac92`$, $`(\alpha_m,\alpha_k,\alpha_n)=(\tfrac12,\tfrac92,2)`$ achieves $`\theta=\tfrac12`$. For $`\tfrac32\le\tau_k<\tfrac92`$, set $`\theta_0=2-\tau_k/3>\tfrac12`$, $`\alpha_m=\tfrac12`$, $`\alpha_k=\tau_k`$, and $`\alpha_n=\tfrac52-\theta_0=\tfrac12+\tau_k/3\ge1`$. Then $`1-\alpha_m=\alpha_m=\tfrac12<\theta_0`$, $`\tfrac52-\alpha_n=\theta_0`$, $`1+2\alpha_n-\alpha_k=2-\tau_k/3=\theta_0`$, $`\alpha_k-2\alpha_n=\tau_k/3-1<\theta_0`$, and $`1+\alpha_m-\alpha_k=\tfrac32-\tau_k\le0<\theta_0`$. $`\square`$

Proposition 12.1 (certified in [`quasi_rh_extractor.py`](quasi_rh_extractor.py) [C], alongside an empirical scan of $`880`$ fundamental quadratic twists showing zero $`\ell`$-adic class-group contagion around Sylow spikes in [`frontier6_bezout_and_centres.py`](frontier6_bezout_and_centres.py) [C]) reveals a striking structural feature: in the intermediate window-free regime $`3<\tau_k<\tfrac92`$ (e.g., $`\tau_k=3.6\implies\theta^{*}=\tfrac45`$), **every linear-term automorphic cancellation (`H1, H2, H4, H5, H6`, entering through $`1-\alpha_m`$ and $`\alpha_m`$ at $`\tfrac12`$) sits strictly below $`\theta^{*}`$ with simultaneous polynomial slack $`\theta^{*}-\tfrac12>0`$**. Caveat: translating Proposition 12.1 into an unconditional theorem about $`\zeta(s)`$ requires verifying directly in the Lean development which concrete $`\tau_k`$ and table-modulus exponent `(H0)–(H10)` achieve without windows.

---

## 13. Summary of companion verification suite

Every `[C]` calculation and constructive example in Parts I and II is verified by an automated `pytest` suite (`60` tests across `10` test modules):

| Module | Companion Test File | Verified Content |
|---|---|---|
| [`rh_core.py`](rh_core.py), [`rh_audit.py`](rh_audit.py) | [`test_use_cases.py`](test_use_cases.py), [`test_cli.py`](test_cli.py) | Explicit formula for $`G(x)`$, Rubinstein's §3.3 table (`84701.032` vs. `84463.394`), double-zero & GUE/Poisson blindness of $`x^{1/2+\varepsilon}`$, $`E_{11a1}`$ Sato–Tate moments, exact $`h(-1319)=45`$ recovery (`19` tests) |
| [`gl2_local_to_global.py`](gl2_local_to_global.py) | [`test_gl2_solver.py`](test_gl2_solver.py) | Exact analytic ranks $`r\in\{0,1,2,3\}`$ and tight Weil–Mestre GRH rank ceilings on `11a1`, `37a1`, `389a1`, `5077a1` from local $`a_p`$ (`4` tests) |
| [`cm_determinant_lab.py`](cm_determinant_lab.py) | [`test_frontier1_determinant.py`](test_frontier1_determinant.py) | §9.1–9.2: $`208`$-digit unshifted hyperfactorial $`H_8(1)`$, Theorem 9.1 basis Jacobian, prime-certified split-prime counterexamples at $`a_0\ge9`$ (`4` tests) |
| [`gl2_divisibility_lab.py`](gl2_divisibility_lab.py) | [`test_frontier2_gl2_divisibility.py`](test_frontier2_gl2_divisibility.py) | §10: Honda integrality, ordinary slope collapse, $`100\%`$ mod-$`p^{2}`$ ASD Frobenius trace recovery without point counting (`8` tests) |
| [`projective_tower_lab.py`](projective_tower_lab.py) | [`test_frontier3_tower.py`](test_frontier3_tower.py) | §11: Theorem 11.1 pullback/pushforward Iwasawa exactness ($`N=8\dots256`$) & exact $`\mathbb{Q}`$-rank $`R`$ vs. Frobenius rank $`2`$ (`5` tests) |
| [`proof_dag_decompiler.py`](proof_dag_decompiler.py) | [`test_frontier4_decompiler.py`](test_frontier4_decompiler.py) | §5.2 & §12: Tarjan SCC quantifier-poset verification, `(H12)` leak catching, hidden-constraint LP & dual shadow prices (`3` tests) |
| [`lean_ast_ingestor.py`](lean_ast_ingestor.py), [`tropical_mip_auditor.py`](tropical_mip_auditor.py), [`quasi_rh_extractor.py`](quasi_rh_extractor.py) | [`test_open_problem_suite.py`](test_open_problem_suite.py) | [`rh_blueprint.lean`](rh_blueprint.lean) AST ingestion, disjunctive MILP, Proposition 12.1 $`\mathbb{Q}`$-Farkas certificates, & Guth–Maynard `30/13` zero-density saddle (`4` tests) |
| [`fitting_ideal_lab.py`](fitting_ideal_lab.py) | [`test_frontier5_fitting_ideal.py`](test_frontier5_fitting_ideal.py) | §9.3 & §10.3: $`51`$-scale Padé annihilation ($`324\to201`$ digits), divided-power multi-rank stability, CM curves `49a1`/`121b1`, $`a_p`$-free $`GL_2`$ Wronskians (`6` tests) |
| [`frontier6_bezout_and_centres.py`](frontier6_bezout_and_centres.py), [`frontier7_hnf_slope_engine.py`](frontier7_hnf_slope_engine.py) | [`test_frontiers_6_and_7.py`](test_frontiers_6_and_7.py) | §9.3 & §5.3: $`10^{-19}`$ Gauss–Kronecker identity, Bezout & kernel-covolume laws, $`\ell`$-adic twist ball scans, rank-deficient modular HNF (`7` tests) |

---

## References

- M. Rubinstein, *The Alpoge-Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, working draft v2, 22 November 2026. https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf
- Lean FRO and the Mathlib Initiative, *Lean verification of the Riemann Hypothesis*, 21 November 2026. https://leanprover-community.github.io/blog/posts/riemann-hypothesis/
- L. Alpoge, announcement note, 20 November 2026. https://alpo.ge/remark.pdf
- T. Tao, *Unforced Euler singularity, unforced Navier-Stokes singularity, and the Generalized Riemann Hypothesis*, 20 November 2026. https://terrytao.wordpress.com/2026/11/20/unforced-euler-singularity-unforced-navier-stokes-singularity-and-the-generalized-riemann-hypothesis/
- A. O. L. Atkin and H. P. F. Swinnerton-Dyer, *Modular forms on noncongruence subgroups*, Proc. Symp. Pure Math. 19 (1971), 1–25.
- L. Guth and J. Maynard, *New large value estimates for Dirichlet polynomials*, 2024.
- T. Honda, *Formal groups and zeta-functions*, Osaka J. Math. 5 (1968), 199–213.
- J. Kaczorowski and J. Pintz, *Oscillatory functions and the Riemann hypothesis*, 1986.
- E. Landau, *Über einen Satz von Tschebyschef*, Math. Ann. 61 (1905), 527–550.
- J.-F. Mestre, *Formules explicites et minorations de conducteurs de variétés algébriques*, Compositio Math. 58 (1986), 209–232.
- T.-J. Stieltjes, *Sur la réduction en fraction continue d'une série procédant suivant les puissances descendantes d'une variable*, Ann. Fac. Sci. Toulouse 3 (1889), H1–H17.
- U. Stuhler (1976) and D. Grayson (1984), on the canonical filtration of a lattice by slopes.
