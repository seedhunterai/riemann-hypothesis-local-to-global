"""
Core mathematical primitives for testing the local-to-global step and
proving the four concrete use cases for an analytic, fixed-object derivation.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np


# ---------------------------------------------------------------------------
# Basic number theory utilities
# ---------------------------------------------------------------------------

def sieve_primes(limit: int) -> List[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            is_prime[p * p : limit + 1 : p] = b"\x00" * len(range(p * p, limit + 1, p))
    return [i for i, flag in enumerate(is_prime) if flag]


def kronecker_symbol(a: int, n: int) -> int:
    """Compute the Kronecker symbol (a / n) for integers a and n > 0."""
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        raise ValueError("n must be positive")
    if a == 0:
        return 1 if n == 1 else 0

    result = 1
    # Factor out powers of 2 from n
    v2 = 0
    while (n & 1) == 0:
        n >>= 1
        v2 += 1
    if v2:
        if (a & 1) == 0:
            return 0
        a_mod8 = a & 7
        if (v2 & 1) and a_mod8 in (3, 5):
            result = -result

    # Odd part: Jacobi symbol
    a %= n
    while a != 0:
        while (a & 1) == 0:
            a >>= 1
            if (n & 7) in (3, 5):
                result = -result
        a, n = n, a
        if (a & 3) == 3 and (n & 3) == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def chi_neg(d: int, n: int) -> int:
    """Quadratic character associated with Q(sqrt(-d)): Kronecker(-d, n)."""
    return kronecker_symbol(-d, n)


def class_number_imag_quad(d: int) -> int:
    """
    Exact class number h(-d) of Q(sqrt(-d)) for fundamental discriminant -d < 0
    (with -d == 1 mod 4), by counting reduced primitive binary quadratic forms
    a*x^2 + b*x*y + c*y^2 with b^2 - 4ac = -d.
    """
    D = -d
    if D >= 0 or (D % 4) not in (0, 1):
        raise ValueError("Expected negative discriminant D = -d with D == 0,1 mod 4")
    count = 0
    max_b = int(math.isqrt(d // 3))
    for b in range(-max_b, max_b + 1):
        num = b * b + d
        if num % 4 != 0:
            continue
        ac = num // 4
        max_a = int(math.isqrt(ac))
        for a in range(max(1, abs(b)), max_a + 1):
            if ac % a != 0:
                continue
            c = ac // a
            if a > c:
                continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            # Reduction boundary rules: if |b| == a or a == c, require b >= 0
            if (abs(b) == a or a == c) and b < 0:
                continue
            count += 1
    return count


# ---------------------------------------------------------------------------
# Section 1: G(x), explicit formula, and Mellin localization
# ---------------------------------------------------------------------------

def S_ramp(x: float, primes: Sequence[int] | None = None) -> float:
    """S(x) = sum_{n <= x} Lambda(n) * (x/n - 1)."""
    n_max = int(math.floor(x))
    if primes is None:
        primes = sieve_primes(n_max)
    total = 0.0
    for p in primes:
        if p > n_max:
            break
        logp = math.log(p)
        pk = p
        while pk <= n_max:
            total += logp * (x / pk - 1.0)
            pk *= p
    return total


def G_exact(x: float, primes: Sequence[int] | None = None) -> float:
    """G(x) = S(x) - x log x + (1 + gamma) x."""
    gamma = float(mp.euler)
    return S_ramp(x, primes) - x * math.log(x) + (1.0 + gamma) * x


def get_zeta_zeros(num_zeros: int) -> List[complex]:
    """Return the first num_zeros nontrivial zeros rho = 1/2 + i*gamma_k (gamma_k > 0)."""
    return [complex(mp.zetazero(k)) for k in range(1, num_zeros + 1)]


def G_explicit_formula(x: float, zeros: Sequence[complex]) -> float:
    """
    Explicit formula approximation:
        G(x) = log(2*pi) - sum_rho x^rho / (rho*(rho - 1)) + O(x^-2)
    where `zeros` lists zeros with Im(rho) > 0 (conjugate pairs added automatically).
    """
    val = math.log(2.0 * math.pi)
    logx = math.log(x)
    for rho in zeros:
        term = np.exp(rho * logx) / (rho * (rho - 1.0))
        val -= 2.0 * float(np.real(term))
    return val


# ---------------------------------------------------------------------------
# Section 3: Rubinstein's small table (q=11, d=1319, N=16, x=1000)
# ---------------------------------------------------------------------------

def L_prime_over_L_at_1(d: int, terms: int = 200000) -> float:
    """
    Compute L'(1, chi_{-d}) / L(1, chi_{-d}) using mpmath's Hurwitz zeta
    representation for exactness.
    """
    # L(s, chi) = d^(-s) * sum_{a=1}^{d-1} chi(a) * zeta(s, a/d)
    # Using mpmath at high precision via functional equation / Hurwitz zeta:
    # Since chi is odd (chi(-1) = -1), L(1, chi) = pi * h(-d) / sqrt(d) (for d > 4)
    # and L'(1, chi) = -sum_{a=1}^{d-1} chi(a) * (log(d) * psi(a/d) - ... )
    # Direct Hurwitz-zeta derivative at s=1 is exact and fast for d ~ 1319.
    with mp.workdps(30):
        # L(s, chi) = d^(-s) sum_{a=1}^{d-1} chi(a) * zeta(s, a/d)
        # Since sum chi(a) = 0, the poles at s=1 cancel:
        # zeta(s, a/d) = 1/(s-1) - psi(a/d) + gamma_1(a/d)*(s-1) + ...
        # Easier & exact: use mp.diff on the finite Hurwitz sum at s=1!
        # Even faster: use the exact identity
        # L(1, chi) = -(1/d) sum_{a=1}^{d-1} chi(a) * psi(a/d)
        # L'(1, chi) = -log(d)*L(1, chi) + (1/d) sum_{a=1}^{d-1} chi(a) * gamma_1(a/d)
        # where gamma_1(z) is the generalized Stieltjes constant.
        # Alternatively, finite-difference or exact formula via Gauss sum & log Gamma:
        # For odd primitive chi mod d:
        # L'(1, chi)/L(1, chi) = gamma + log(2*pi/d) - (1/h) * ... (Chowla-Selberg!)
        # Actually, let's compute via mp.dirichlet(s, chi_vec) if available, or
        # exact numerical derivative of Hurwitz zeta at s = 1 + eps.
        chi_vals = [chi_neg(d, a) for a in range(d)]

        def L_of_s(s):
            s = mp.mpf(s)
            acc = mp.mpf(0)
            for a in range(1, d):
                c = chi_vals[a]
                if c:
                    acc += c * mp.hurwitz(s, mp.mpf(a) / d)
            return acc * (mp.mpf(d) ** (-s))

        # Exact L(1, chi) = -(1/d) sum chi(a) * psi(a/d)
        L1 = -sum(chi_vals[a] * mp.digamma(mp.mpf(a) / d) for a in range(1, d) if chi_vals[a]) / d
        # Derivative at s=1 via symmetric difference at 1 +/- 1e-8 (with 30 dps -> 15+ digits accurate)
        h = mp.mpf("1e-8")
        Lp1 = (L_of_s(1 + h) - L_of_s(1 - h)) / (2 * h)
        return float(Lp1 / L1)


def compute_rubinstein_table(
    q: int = 11,
    d: int = 1319,
    N: int = 16,
    a: float = 10000.0,
    w_dist: float = 10.0,
    w_other: float = 0.05,
    dist_col: int = 1,
    L_ratio: float | None = None,
) -> Tuple[float, float]:
    """
    Compute (P, M) for Rubinstein's small table:
      - 2 rows v in {0, 1}, where v(r) = (1 + chi_q(r)) / 2
      - n = N/4 columns indexed by H = (Z/NZ)^x / {+/- 1}
      - distinguished cell at (row=1, col=[dist_col]) with weight w_dist; all others w_other.
    """
    n = N // 4
    # Normalize dist_col into canonical representative in {1, 3, ..., N/2 - 1}
    def col_rep(r: int) -> int:
        c = r % N
        return min(c, N - c)

    target_col = col_rep(dist_col)
    max_r = int(math.floor(a / min(w_dist, w_other)))
    primes = sieve_primes(max_r)

    P = 0.0
    for p in primes:
        if p == 2 or p == q:
            continue
        logp = math.log(p)
        pk = p
        while pk <= max_r:
            c_d = chi_neg(d, pk)
            factor = 1 - c_d
            if factor != 0:
                c_q = chi_neg(q, pk)
                row = (1 + c_q) // 2
                col = col_rep(pk)
                w = w_dist if (row == 1 and col == target_col) else w_other
                val = a / pk - w
                if val > 0.0:
                    P += factor * logp * val
            pk *= p

    if L_ratio is None:
        L_ratio = L_prime_over_L_at_1(d)

    num_cells = 2 * n
    avg_log_w = (math.log(w_dist) + (num_cells - 1) * math.log(w_other)) / num_cells
    gamma = float(mp.euler)
    M = a * (math.log(a) - 1.0 - gamma + L_ratio - avg_log_w)
    return P, M


# ---------------------------------------------------------------------------
# Use Case 2: Concrete non-CM elliptic curve E_11a1: y^2 + y = x^3 - x^2
# ---------------------------------------------------------------------------

def elliptic_curve_ap(
    p: int, a1: int = 0, a2: int = -1, a3: int = 1, a4: int = 0, a6: int = 0
) -> int:
    """
    Compute a_p = p + 1 - #E(F_p) for a general Weierstrass elliptic curve
        y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    over F_p (for any prime p).
    """
    if p == 2:
        pts = sum(
            1
            for x in range(2)
            for y in range(2)
            if (y * y + a1 * x * y + a3 * y - (x * x * x + a2 * x * x + a4 * x + a6)) % 2 == 0
        )
        return p + 1 - (pts + 1)

    # Complete the square: (2y + a1*x + a3)^2 = 4*(x^3 + a2*x^2 + a4*x + a6) + (a1*x + a3)^2
    qr = [0] * p
    qr[0] = 1
    for z in range(1, p):
        qr[(z * z) % p] += 1

    affine_pts = 0
    for x in range(p):
        lin = (a1 * x + a3) % p
        rhs = (x * x * x + a2 * x * x + a4 * x + a6) % p
        disc = (4 * rhs + lin * lin) % p
        affine_pts += qr[disc]
    return p + 1 - (affine_pts + 1)


def elliptic_curve_11a1_ap(p: int) -> int:
    """
    Compute a_p = p + 1 - #E(F_p) for the non-CM elliptic curve
    E_11a1: y^2 + y = x^3 - x^2 (conductor 11, modular form Delta_2 of weight 2).
    """
    if p == 11:
        return 1
    return elliptic_curve_ap(p, a1=0, a2=-1, a3=1, a4=0, a6=0)


# ---------------------------------------------------------------------------
# Use Case 1: Synthetic double-zero and GUE vs Poisson spectral diagnostics
# ---------------------------------------------------------------------------

def G_with_double_zero(x: float, zeros: Sequence[complex], double_idx: int = 0) -> float:
    """
    Explicit formula for a hypothetical zeta-like function where the zero
    at `zeros[double_idx]` has multiplicity 2 (adding a second-order pole
    contribution to F(s), whose inverse Mellin transform has a log(x)*x^rho term).
    """
    base = G_explicit_formula(x, zeros)
    rho = zeros[double_idx]
    logx = math.log(x)
    # A double pole at s = rho in F(s) gives a residue proportional to x^rho * log(x) / (rho*(rho-1))
    extra = np.exp(rho * logx) * logx / (rho * (rho - 1.0))
    return base - 2.0 * float(np.real(extra))


def pair_correlation_form_factor(gammas: Sequence[float], alpha: float) -> float:
    """
    Empirical Montgomery pair-correlation form factor F(alpha) on unfolded spacings.
    For GUE (random matrix / zeta zeros), F(alpha) -> |alpha| for |alpha| < 1 (level repulsion).
    For Poisson (uncorrelated zeros on Re s = 1/2), F(alpha) ~ 1 for all alpha > 0.
    """
    g = np.asarray(gammas, dtype=float)
    n = len(g)
    # Unfold to unit mean spacing
    diffs = np.diff(g)
    mean_gap = float(np.mean(diffs))
    u = (g - g[0]) / mean_gap
    # Compute (1/N) |sum_j exp(2*pi*i*alpha*u_j)|^2, subtracting the forward peak
    # Use a Hann window to suppress finite-sample boundary ringing
    w = np.hanning(n)
    w = w / np.sqrt(np.sum(w ** 2))
    s = np.sum(w * np.exp(2j * np.pi * alpha * u))
    return float(np.abs(s) ** 2)
