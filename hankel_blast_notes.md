# Beyond BLAST: Hankel-Transformed Windows

## Original BLAST Structure

In standard BLAST the angular power spectrum is:

\[
C_\ell^{AB} = \int dk\; k^2 \int d\chi_1 \int d\chi_2\; W^A(\chi_1)\; W^B(\chi_2)\; j_\ell(k\chi_1)\; j_\ell(k\chi_2)\; P(k)
\]

The key observation is that the \(\chi\) integrals are **easy** — the windows \(W^A(\chi)\) are smooth, slowly varying functions, so simple quadrature suffices. The \(k\) integral is **hard** because of the two rapidly oscillating Bessel functions \(j_\ell(k\chi_1)\, j_\ell(k\chi_2)\).

BLAST solves this by:

1. Chebyshev-expanding \(P(k)\) → coefficients \(p_n\)
2. Precomputing
\[
\tilde{T}_\ell^{AB}(\chi_1, \chi_2) = \int dk\; k^2\; T_n(k)\; j_\ell(k\chi_1)\; j_\ell(k\chi_2)
\]
once and for all
3. Contracting:
\[
C_\ell^{AB} = \sum_n p_n \int d\chi_1\, d\chi_2\; W^A(\chi_1)\, W^B(\chi_2)\; \tilde{T}_n^{AB}(\chi_1, \chi_2)
\]

**Split:** hard \(k\)-integral precomputed → easy \(\chi\)-integrals done at runtime.

---

## New Structure with Hankel-Transformed Windows

The angular power spectrum in SFB space is:

\[
S_\ell^{AB}(k_1, k_2) = \int dk\; k^2\; P(k)\; \widetilde{W}_\ell^A(k, k_1)\; \widetilde{W}_\ell^B(k, k_2)
\]

where each Hankel-transformed window is itself an integral:

\[
\widetilde{W}_\ell^A(k, k_1) = \int d\chi\; \chi^2\; f^A(\chi)\; j_\ell(k\chi)\; j_\ell(k_1\chi)
\]

Expanding everything, the full object is a **triple integral**:

\[
S_\ell^{AB}(k_1, k_2) = \int dk\; k^2\, P(k)
\int d\chi\; \chi^2\, f^A(\chi)\; j_\ell(k\chi)\; j_\ell(k_1\chi)
\int d\chi'\; \chi'^2\, f^B(\chi')\; j_\ell(k\chi')\; j_\ell(k_2\chi')
\]

---

## The Critical Inversion Relative to BLAST

| | Original BLAST | New code |
|---|---|---|
| \(\chi\) integrals | **Easy** — smooth windows, no Bessel | **Hard** — two Bessel functions per window |
| \(k\) integral | **Hard** — two Bessels, needed Chebyshev | **Easy** — just \(k^2 P(k)\), no Bessel |

---

## What You Precompute and What You Contract

The Chebyshev strategy now applies to the **\(\chi\) integrals**, not the \(k\) integral. For each probe, precompute:

\[
\widetilde{W}_\ell^A(k, k_1) = \int d\chi\; \chi^2\; f^A(\chi)\; j_\ell(k\chi)\; j_\ell(k_1\chi)
\]

This table has shape `(n_k, n_k1)` — both \(k\) and \(k_1\) live on the Chebyshev grid. It is computed **once per probe per \(\ell\)** and stored.

The final \(k\)-convolution is then trivial:

\[
S_\ell^{AB}(k_1, k_2) = \int dk\; k^2\; P(k)\; \widetilde{W}_\ell^A(k, k_1)\; \widetilde{W}_\ell^B(k, k_2)
\]

For fixed \((k_1, k_2)\) this is a smooth 1D integral — one slice from each table, multiplied by \(k^2 P(k)\), integrated with standard quadrature. No Chebyshev decomposition of \(P(k)\) is needed here.

---

## Analogy with `compute_T̃`

The original `compute_T̃` computes, for each \((\chi, R)\) pair:

\[
\tilde{T}[\ell, \chi, R, n] = \int dk\; k^\beta\; T_n(\log k)\; j_\ell(k\chi)\; j_\ell(kR\chi)
\]

where:

- Integration variable: \(k\), on a Clenshaw-Curtis grid of \(N = 2^{15}+1\) points
- `Bessel1[i,k]` \(= j_\ell(k \cdot \chi_i)\)
- `Bessel2[i,k]` \(= j_\ell(k \cdot R \cdot \chi_i)\)
- `T[n,k]` = \(n\)-th Chebyshev polynomial evaluated at \(\log k\)
- Output shape: `(1, nχ, nR, n_cheb+1)`

The power spectrum is reconstructed as \(P(k) = \sum_n c_n\, T_n(\log k)\), so the contraction `w_ell_tullio(c, T_tilde)` gives:

\[
w[\ell, \chi, R] = \int dk\; k^\beta\, P(k)\, j_\ell(k\chi)\, j_\ell(kR\chi)
\]

---

## Step 1 — The New Integral to Precompute

For each probe \(A\), each tomographic bin \(b\), each \(\ell\):

\[
\widetilde{W}_\ell^{A,b}(k, k') = \int d\chi\; \chi^2\; f^{A,b}(\chi)\; j_\ell(k\chi)\; j_\ell(k'\chi)
\]

- Integration variable: \(\chi\) (not \(k\))
- Both \(k\) and \(k'\) live on the **same Chebyshev \(k\)-grid** used in the original code
- Output: a 2D matrix indexed by \((k, k')\), shape `(n_k, n_k)` per bin per \(\ell\)
