# Roadmap for Implementing $S_l^{den}(k_1, k_2)$

This roadmap is a working guide for implementing the density correlation equation in `equations.tex` and then wiring it into `beyond_BLAST.ipynb`.

Reference paper: BLAST, arXiv:2410.03632v1, especially the discussion of the 3x2pt algorithm, the Chebyshev factorization of the power spectrum, and the N5K validation flow.

## 1. Target equation

The quantity to implement is

$$
S_l^{den}(k_1,k_2) = N_l^{den} \int dk\, k^2 P_{lin}(k)
\widetilde{W}_l^{den}(k_1,k)\widetilde{W}_l^{den}(k_2,k),
$$

with

$$
\widetilde{W}_l^{den}(k_1,k) = \int d\chi\, \chi^2 \, \frac{H(\chi)b(\chi)n(\chi)D(\chi)}{c}\, j_l(k\chi)j_l(k_1\chi).
$$

The key idea is to keep the structure consistent with the BLAST notebook: first build background quantities, then build kernels, then perform the final wavenumber integral.

## 2. Reuse what already exists

Before adding any new code, reuse the pieces that are already present in `blast_code/src` and shown in `BLAST_x_N5K.ipynb`.

- `compute_adimensional_hubble_factor(z, cosmo)` in `blast_code/src/background.jl` already gives the normalized Hubble factor `E(z)`.
- `compute_χ(z, cosmo)` in `blast_code/src/background.jl` already provides comoving distance.
- `compute_hubble_factor(z, cosmo)` in `blast_code/src/background.jl` already gives `H(z)`.
- `evaluate_background_quantities!` in `blast_code/src/background.jl` already fills the background arrays used later by kernels.
- `compute_kernel!` in `blast_code/src/background.jl` already constructs the N5K-style galaxy and shear kernels.
- `compute_T̃`, `fast_chebcoefs`, `w_ell_tullio`, and `compute_Cℓ` in `BLAST_x_N5K.ipynb` show the current BLAST workflow from power-spectrum decomposition to final angular power spectra.

If a function already exists in the source tree, use it directly instead of introducing a duplicate name in the notebook.

## 3. What needs to be added

The density projection in the roadmap requires three new conceptual pieces.

1. A density weight function
   - Define `f^{den}(\chi) = H(\chi)b(\chi)n(\chi)D(\chi)/c`.
   - If the notebook already has `b(z)`, `n(z)`, and `D(z)` as arrays or interpolants, wrap them in a small callable function.

2. A density window function
   - Define `W_l^{den}(k,\chi) = f^{den}(\chi) j_l(k\chi)`.
   - This is the direct analogue of the BLAST window-building stage.

3. A Hankel transform and final `S_l^{den}` integration
   - Compute `\widetilde{W}_l^{den}(k_1,k)` by integrating over `\chi`.
   - Then integrate over `k` with `P_{lin}(k)` to obtain `S_l^{den}(k_1,k_2)`.

## 4. Suggested implementation path in `beyond_BLAST.ipynb`

### Step 1: make the background explicit

Use the existing background functions from `blast_code/src/background.jl` to evaluate:

- `E(z)`
- `H(z)`
- `\chi(z)`

This is already consistent with the notebook setup in `BLAST_x_N5K.ipynb`, where `z_b` and `\chi_b` are loaded and interpolated.

### Step 2: define the density weights

Create a small notebook cell that defines `f_den(chi)` or `f_den(z)` using the same conventions as the paper and the BLAST notebook.

Recommended ingredients:

- `H(z)` from the source tree
- galaxy bias `b(z)`
- selection function `n(z)`
- growth factor `D(z)` if you are using the full density expression from the equations file

### Step 3: implement the density window function

Add a function that mirrors the notebook logic used for the clustering and lensing kernels:

```julia
W_den(l, k, chi) = f_den(chi) * spherical_bessel(l, k * chi)
```

If `spherical_bessel` is not yet available in `blast_code/src`, implement it first or inline it locally for the notebook prototype.

### Step 4: implement `\widetilde{W}_l^{den}`

Use numerical quadrature over `\chi`.

Practical notes:

- reuse the same `\chi` grid that is already used in the notebook
- keep the integration range consistent with the BLAST notebook and the N5K setup
- cache the result if you evaluate many `(k_1, k)` pairs

### Step 5: implement `S_l^{den}(k_1,k_2)`

Follow the same pattern as the BLAST notebook's final `C_l` construction:

1. build a `k` grid
2. evaluate `P_lin(k)`
3. evaluate `\widetilde{W}_l^{den}(k_1,k)` and `\widetilde{W}_l^{den}(k_2,k)`
4. integrate the product over `k`

## 5. Validation checks

The notebook should verify the implementation in small, cheap steps before any expensive grid run.

1. Background check
   - confirm that `E(z)` and `\chi(z)` match the precomputed arrays already used in `BLAST_x_N5K.ipynb`

2. Kernel sanity check
   - plot `f^{den}(\chi)` and `W_l^{den}(k,\chi)` for a few representative values of `l` and `k`

3. Transform check
   - inspect `\widetilde{W}_l^{den}(k_1,k)` for smoothness and expected oscillatory structure

4. Limber comparison
   - compare the final result against the Limber-limit construction already used in `BLAST_x_N5K.ipynb`

5. Scaling check
   - measure runtime as the `k` and `\chi` grids grow, so the notebook remains usable for larger parameter scans

## 6. Notebook order of operations

The best order in `beyond_BLAST.ipynb` is:

1. load or compute background quantities
2. define density weights
3. define spherical Bessel helpers if needed
4. compute `W_l^{den}`
5. compute `\widetilde{W}_l^{den}`
6. compute `S_l^{den}(k_1,k_2)`
7. compare with the Limber-style approximation used in the existing notebook

## 7. Where to edit

- `sheet/equations.tex` for the derivation and final displayed equations
- `beyond_BLAST.ipynb` for the notebook implementation
- `blast_code/src/background.jl` only if a missing helper truly belongs in the shared source tree

## 8. Expected outcome

After this roadmap is followed, the notebook should contain a reproducible path from the background model to `S_l^{den}(k_1,k_2)`, with the implementation aligned to the BLAST paper and with the source tree reused whenever an equivalent function already exists.