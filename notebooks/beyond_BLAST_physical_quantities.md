# Physical Quantities Reference: `beyond_BLAST` Notebook and Supporting Functions

---

## Overview

This notebook computes angular power spectra \( C_\ell \) for 3×2pt cosmological observables — galaxy–galaxy, shear–shear, and galaxy–shear cross-correlations — using a non-Limber approach based on Chebyshev polynomial decomposition (the BLAST/beyond-BLAST framework). The computation proceeds from background cosmology through observable kernels, spectral decomposition, and final power spectrum assembly.

---

## 1. Background Cosmology and Coordinate Variables

### 1.1 Redshift \( z \) and Comoving Distance \( \chi \)

- `z_b`, `z_range`, `zed`: Arrays of **cosmological redshift** \( z \), a dimensionless quantity measuring how much the universe has expanded since light was emitted. \( z = 0 \) is today; higher values correspond to earlier cosmic times.
- `χ_b`, `x_range`, `chi`: Arrays of **comoving distance** \( \chi \), measured in \(\mathrm{Mpc}/h\). The comoving distance is defined as:
  \[
  \chi(z) = \int_0^z \frac{c \, dz'}{H(z')}
  \]
  where \( c \) is the speed of light and \( H(z') \) is the Hubble parameter. It represents the distance to a source factoring out the expansion of the universe, so that comoving positions of objects that move only with the Hubble flow stay fixed. In the code, `x_range` spans from 26 to 7000 Mpc/h (96 points), and `x_min`, `x_max` are its bounds.
- `z_of_χ`, `chi_of_z`: **Akima interpolation** objects mapping \( \chi \mapsto z(\chi) \) and \( z \mapsto \chi(z) \) respectively. The Akima interpolation is a local cubic spline that is robust to outliers and does not exhibit the oscillatory behavior of global polynomial interpolants. It is used throughout to convert between the two coordinate variables.
- `zmin`, `zmax`: Minimum and maximum redshift of the galaxy redshift distribution, loaded from the N5K data file (`dNdzs_fullwidth.npz`).
- `n_chi` = 96: Number of discrete comoving distance (or redshift) grid points used to represent radial kernels.

### 1.2 Cosmological Model

- `cosmo`: A **flat ΛCDM cosmology** object (`Blast.FlatΛCDM()`), containing the standard parameters: Hubble constant \( H_0 \), matter density parameter \( \Omega_m \), dark energy density \( \Omega_\Lambda \), etc.
- `C_LIGHT`: Speed of light \( c \) in km/s (numerical value \( \approx 2.998 \times 10^5 \) km/s).

---

## 2. The Hubble Function and Growth Factor

### 2.1 Dimensionless Hubble Factor \( E(z) \)

Computed via `Blast.compute_adimensional_hubble_factor(z, cosmo)`. This is:
\[
E(z) \equiv \frac{H(z)}{H_0} = \sqrt{\Omega_m (1+z)^3 + \Omega_\Lambda}
\]
for a flat ΛCDM cosmology. It is dimensionless and encodes how the expansion rate changes with redshift relative to its present-day value.

### 2.2 Physical Hubble Factor \( H(z) \)

Computed via `Blast.compute_hubble_factor(z, cosmo)`, which returns \( H(z) = H_0 \cdot E(z) \) in units of \(\mathrm{km\,s^{-1}\,Mpc^{-1}}\). In the galaxy kernel, the ratio \( c/H(z) \) — stored in `inv_Hubble_array` — serves as a Jacobian factor converting from redshift to comoving distance:
\[
\frac{c}{H(z)} = \frac{d\chi}{dz}
\]
with units of Mpc/h. This factor appears because the radial integrals are performed over redshift rather than comoving distance.

### 2.3 Heath Integral and Linear Growth Factor \( D(z) \)

The **linear growth factor** \( D(z) \) describes how small density perturbations grow over time in the linear regime. The code uses the Heath (1977) integral formula:
\[
D^{\mathrm{unnorm}}(z) = E(z) \int_z^{\infty} \frac{1+z'}{E(z')^3} \, dz'
\]
implemented in the function `heath_integral(cosmo, z)`, which calls `quadgk` (Gauss–Kronrod adaptive quadrature) from \( z \) to \( \infty \) with relative tolerance \( 10^{-10} \).

The **normalized growth factor** is:
\[
D(z) = \frac{D^{\mathrm{unnorm}}(z)}{D^{\mathrm{unnorm}}(0)}
\]
so that \( D(0) = 1 \). It is stored in `D_growth_array` and `growth`. This normalization ensures that the amplitude of the power spectrum at \( z = 0 \) is taken as the reference, and the factor \( D(z) \) tells you by how much the fluctuation amplitude was smaller at redshift \( z \).

---

## 3. Source Redshift Distribution \( n(z) \)

The **normalized source redshift distribution** \( n(z) \) models the probability distribution of observed galaxy redshifts. Its un-normalized form is:
\[
n(z) = A \left(\frac{z}{z_0}\right)^\alpha \exp\!\left[-\left(\frac{z}{z_0}\right)^\beta\right]
\]
with parameters fixed in the code as:
- \( z_0 = 0.9/\sqrt{2} \approx 0.636 \)
- \( A = 1.5/z_0 \)
- \( \alpha = 2 \)
- \( \beta = 1.5 \)

This is the standard Smail-type parametrization used for LSST-like surveys. The normalized version `nz_norm` satisfies:
\[
\int n(z) \, dz = 1
\]
achieved by dividing by the integral computed with Akima interpolation and `quadgk`. The normalization makes \( n(z) \) a proper probability density over redshift, so it can be physically interpreted as the fraction of sources per unit redshift interval.

---

## 4. Galaxy Bias \( b(z) \)

The **linear galaxy bias** \( b(z) \) relates the galaxy overdensity field to the underlying dark matter overdensity:
\[
\delta_g(\mathbf{x}, z) = b(z) \, \delta_m(\mathbf{x}, z)
\]
The code uses:
\[
b(z) = b_0 \sqrt{1+z}, \quad b_0 = 1.0
\]
stored in `bz_array` (and passed as `bias`). This simple model reflects the observed mild redshift evolution of galaxy bias, increasing with \( z \) because galaxies at higher redshift tend to inhabit more massive, rarer haloes.

---

## 5. Window (Kernel) Functions

The window functions encode how each observable probe is weighted along the line of sight. They directly enter the projection integrals that produce the angular power spectra.

### 5.1 Galaxy–Galaxy Kernel \( W^g(z) \)

Computed in `galaxy_prefactor` and stored in `prefac` / `gal_prefact_W`:
\[
W^g(z) = \chi(z)^2 \, b(z) \, D(z) \, n(z)
\]
The factor \( \chi(z)^2 \) is the **comoving volume element** per steradian (i.e., the Jacobian that turns a redshift-space integral into a solid-angle weighted comoving-volume integral). Together, the four factors \( \chi^2 b D n \) weight the contribution of shells at redshift \( z \) to the observed galaxy angular clustering: the more sources there are, the stronger the bias, and the larger the growth factor, the stronger the signal at that redshift. Units are approximately \((\mathrm{Mpc}/h)^2\).

### 5.2 Lensing Efficiency \( q(z) \)

Computed in `compute_lensing_efficiency` and stored in `lens_int_array`:
\[
q(z_l) = \int_{z_l}^{z_{\max}} n(z_s) \left(1 - \frac{\chi(z_l)}{\chi(z_s)}\right) dz_s
\]
This integral is computed with `quadgk` (relative tolerance \( 10^{-5} \)) for each lens redshift \( z_l \). The integrand is zero whenever the source distance \( \chi(z_s) \leq \chi(z_l) \) (a source cannot be lensed by a structure at greater distance than itself). The factor \( (1 - \chi_l/\chi_s) \) is the **lensing geometric efficiency**: it is zero when lens and source are at the same distance and approaches 1 for sources at infinity. The integral therefore accumulates contributions from all sources behind the lens, weighted by their geometric lensing efficiency.

### 5.3 Shear–Shear Kernel \( W^\kappa(z) \)

Computed in `shear_prefactor` and stored in `prefac_shear` / `shear_prefac_W`:
\[
W^\kappa(z) = \frac{1}{\chi(z)^2} \cdot \frac{c}{H(z)} \cdot D(z) \cdot \frac{3 H_0^2 \Omega_m}{2 c^2} \cdot (1+z) \cdot q(z)
\]
where `pref = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2` = \( \frac{3 H_0^2 \Omega_m}{2c^2} \) is the standard **lensing prefactor** in units of \( (\mathrm{Mpc}/h)^{-2} \). This factor comes from the Poisson equation relating the gravitational potential to the matter overdensity and from the definition of the convergence \( \kappa \). The factor \( (1+z) \) converts from comoving to physical distances. The factor \( \chi^{-2} \) comes from the angular diameter distance in the lensing Jacobian. The factor \( c/H(z) \) is again the comoving distance–redshift Jacobian.

---

## 6. Chebyshev Representation of Kernels

### 6.1 Clenshaw–Curtis Grid and Weights

The code maps all radial integrals onto a **Clenshaw–Curtis (CC) quadrature** grid. This is an optimal quadrature rule based on the zeros of Chebyshev polynomials, chosen because it achieves spectral convergence (exponential error reduction with increasing node count) for smooth functions.

- `get_clencurt_grid_z(zmin, zmax, N)`: Returns the \( N \) Clenshaw–Curtis nodes mapped from the standard interval \([-1, 1]\) to \([z_{\min}, z_{\max}]\) via the affine transformation:
  \[
  z_m = \frac{z_{\max}-z_{\min}}{2} \cos\left(\frac{m\pi}{N-1}\right) + \frac{z_{\min}+z_{\max}}{2}
  \]
  The endpoints are nudged by \( \pm 10^{-8} \) to avoid exact boundary issues. Stored in `z_cheb_nodes`.

- `get_clencurt_weights_z(zmin, zmax, N)`: Returns the corresponding integration weights, scaled from \([-1,1]\) to \([z_{\min}, z_{\max}]\) as:
  \[
  w_m^{(z)} = \frac{z_{\max}-z_{\min}}{2} w_m^{CC}
  \]
  where \( w_m^{CC} \) are the standard Clenshaw–Curtis weights. These weights implement the quadrature rule:
  \[
  \int_{z_{\min}}^{z_{\max}} f(z) \, dz \approx \sum_{m=1}^{N} w_m^{(z)} f(z_m)
  \]
  Stored in `w` (local), `w_k`, `w_kp`, `w_kpp` (for \( k \)-space grids).

- `n_cheb` = 119: **Number of Chebyshev nodes** (polynomial degree + 1), i.e., the number of Chebyshev basis functions used to expand the kernel.
- `N` = \( 2^{15}+1 \): **Number of CC quadrature points** in \( k \)-space used to evaluate the Bessel-function integrals. This large value is needed because spherical Bessel functions are highly oscillatory.

### 6.2 Chebyshev Polynomial Matrix \( T_{n,m} \)

Computed in `bessel_cheb_eval_beyond`. The matrix `T` has shape `(n_cheb+1, N)`, where entry \( T_{n,m} \) is the value of the \( n \)-th Chebyshev polynomial of the first kind at the \( m \)-th CC node (mapped to \([-1,1]\)):
\[
x_m = \frac{2 z_m - (z_{\max}+z_{\min})}{z_{\max}-z_{\min}} \in [-1, 1]
\]
\[
T_0(x) = 1, \quad T_1(x) = x, \quad T_n(x) = 2x \, T_{n-1}(x) - T_{n-2}(x)
\]
The three-term recurrence is used for efficiency, avoiding calling trigonometric functions repeatedly. The Chebyshev polynomials form an orthogonal basis on \([-1,1]\) under a particular weight, and any smooth function can be expanded in this basis with rapidly decaying coefficients.

### 6.3 Spherical Bessel Function Matrix

Computed in `bessel_cheb_eval_beyond`. The matrix `Bessel` has shape `(Nz, N)`, where entry \( j_\ell(\chi_i k_m) \) is the **spherical Bessel function of the first kind of order \( \ell \)**, evaluated at the product of the \( i \)-th comoving distance grid point and the \( m \)-th \( k \)-mode. Spherical Bessel functions appear because the plane-wave basis \( e^{i\mathbf{k}\cdot\mathbf{r}} \) decomposes into partial waves on a sphere as:
\[
e^{i\mathbf{k}\cdot\mathbf{r}} = 4\pi \sum_{\ell m} i^\ell j_\ell(kr) Y_{\ell m}^*(\hat{k}) Y_{\ell m}(\hat{r})
\]
The computation is parallelized over \( k \)-columns to exploit Julia's column-major memory layout.

### 6.4 Kernel Values on Chebyshev Nodes

`W_vals` in `galaxy_prefactor_cheb` and `shear_prefactor_cheb`: The window function \( W(z) \) evaluated on the Clenshaw–Curtis nodes `z_cheb_nodes`. For the galaxy case:
\[
W_{\text{vals},m} = \chi(z_m)^2 \, b(z_m) \, D(z_m) \, n(z_m)
\]
All quantities are obtained by Akima interpolation from the original `z_range` to the Chebyshev nodes.

### 6.5 Chebyshev Coefficients

`cheb_coeff_gal`, `cheb_coeff_shear`: The **Chebyshev expansion coefficients** of the window functions, computed via a Fast Fourier Transform (`Blast.fast_chebcoefs`). If \( W(z) \approx \sum_{n=0}^{N_c} c_n T_n(z) \), the array `cheb_coeff` stores the \( c_n \). This allows the smooth part of the integrand (the kernel) to be analytically represented, leaving only the highly oscillatory Bessel functions to be handled numerically.

---

## 7. The \( \tilde{W} \) Tensor and the \( W_{\text{final}} \) Tensor

### 7.1 \( \tilde{W}^{(\ell)}_{i,p,n} \): The Chebyshev–Bessel Overlap Matrix

Computed in `compute_W_tilde`, stored in `W_tilde` with shape `(N_\ell, N_\chi, N_\chi, N_\text{cheb})` = `(100, 96, 96, 119)`. The mathematical definition:
\[
\tilde{W}^{(\ell)}_{i,p,n} = \sum_{m=1}^{N_k} w_m^{(k)} \, T_n(k_m) \, j_\ell(\chi_i k_m) \, j_\ell(\chi_p k_m)
\]
This quantity is the central object of the beyond-BLAST algorithm. It tells you, for a given multipole \( \ell \), how strongly two radial shells at comoving distances \( \chi_i \) and \( \chi_p \) are coupled through the \( n \)-th Chebyshev mode of the kernel. The sum over \( k \) is a Clenshaw–Curtis quadrature approximation to the integral:
\[
\tilde{W}^{(\ell)}_{i,p,n} \approx \int_{k_{\min}}^{k_{\max}} T_n(k) \, j_\ell(\chi_i k) \, j_\ell(\chi_p k) \, dk
\]

In the optimized implementation, since \( \tilde{W} \) does not depend on \( p \) (both Bessel arguments use the same grid), the matrix multiplication `C = A * T'` is used, where `A[i, m] = Bessel[i,m]^2 * w[m]` has shape `(Nk, N)` and `T'` has shape `(N, n_cheb+1)`. The result `C` has shape `(Nk, n_cheb+1)`. The full tensor is filled by tiling `C` along the `p` dimension.

### 7.2 \( W_{\text{final}} \): The Contracted Kernel Tensor

Computed by the `@tullio` contraction:
```julia
W_final_gal[i, j, k] := W_tilde[i, j, k, l] * cheb_coeff_gal[l]
W_final_shear[i, j, k] := W_tilde[i, j, k, l] * cheb_coeff_shear[l]
```
Shape: `(N_\ell, N_\chi, N_\chi)` = `(100, 96, 96)`.

This contraction implements:
\[
W_{\text{final}}^{(\ell)}_{i,p} = \sum_{n=0}^{N_c} c_n \, \tilde{W}^{(\ell)}_{i,p,n} = \sum_{n=0}^{N_c} c_n \sum_{m} w_m T_n(k_m) j_\ell(\chi_i k_m) j_\ell(\chi_p k_m)
\]
Substituting the Chebyshev expansion \( W(k) \approx \sum_n c_n T_n(k) \), this becomes:
\[
W_{\text{final}}^{(\ell)}_{i,p} \approx \int_{k_{\min}}^{k_{\max}} W(k) \, j_\ell(\chi_i k) \, j_\ell(\chi_p k) \, dk
\]
i.e., \( W_{\text{final}} \) is the \( k \)-integrated Bessel transform of the kernel \( W(k) \), evaluated at all pairs \( (\chi_i, \chi_p) \).

---

## 8. The Matter Power Spectrum \( P(k, z) \)

### 8.1 Data and Interpolation

- `Pklin`, `k`, `z`: Arrays of **linear matter power spectrum** values \( P^{\rm lin}(k, z) \), wavenumbers \( k \) (in \( h/\mathrm{Mpc} \)), and redshifts.
- `Pknonlin`: **Non-linear matter power spectrum** \( P^{\rm nl}(k, z) \), typically computed with a fitting formula (e.g., HaloFit).

Both are interpolated in 2D using bicubic B-spline interpolation (`Interpolations.jl`) in \(\log_{10} k\) and linear \( z \), operating in log-space since \( P(k) \) spans many orders of magnitude:
```julia
InterpPmm   # for P^lin
InterpPmm_nl  # for P^nl
```
Extrapolation uses linear extension outside the tabulated range.

### 8.2 Power Spectrum Callables

```julia
power_spectrum(k, χ1, χ2)    = sqrt(P(z(χ1), k) * P(z(χ2), k))
power_spectrum_nl(k, χ1, χ2) = sqrt(P_nl(z(χ1), k) * P_nl(z(χ2), k))
```
These return \( \sqrt{P(k,z(\chi_1)) P(k,z(\chi_2))} \). This geometric-mean form appears in the non-Limber decomposition of the power spectrum as a separable approximation to \( P(k, z_{\rm eff}) \), allowing the covariance between the two radial shells to be factored symmetrically. In practice `power_spectrum_nl` is used in the final computation.

---

## 9. Wavenumber Grid and Integration Weights

- `kmin` = \( 2.5/\chi_{\max} \) in \( h/\mathrm{Mpc} \): Sets the **minimum wavenumber**, corresponding to modes with wavelength comparable to the survey depth — modes larger than the survey volume cannot be meaningfully constrained.
- `kmax` = \( 200/13 \approx 15.4 \, h/\mathrm{Mpc} \): Sets the **maximum wavenumber**, corresponding to the smallest scales included in the analysis. This cutoff removes modes where non-linear effects are too strong to be described by the perturbation-theory-calibrated power spectrum.
- `k_grid`, `kp_grid`, `kpp_grid`: Three independent Clenshaw–Curtis wavenumber grids of length `Nz = n_chi = 96`, each spanning \([k_{\min}, k_{\max}]\). Three grids are defined because the final power spectrum involves a triple \( k \)-integral (over \( k \), \( k' \), \( k'' \)), and having named copies makes the contraction structure transparent.
- `w_k`, `w_kp`, `w_kpp`: Corresponding Clenshaw–Curtis quadrature weights for the three grids.
- `k_cheb`: Chebyshev nodes in \(\log_{10} k\) space (not in \( k \) directly), used to represent the power spectrum as a Chebyshev series in log scale.

---

## 10. Integration Weights with Power Spectrum: `weight_gal`, `weight_shear`, `weight_gal_shear`

These three arrays absorb the \( k \)-dependent part of the power spectrum integration into the quadrature weights. They are computed as:

```julia
Pk_grid         = power_spectrum_nl.(k_grid, first(x), last(x))
weight_gal      = w_k.^2 .* k_grid.^2 .* Pk_grid
weight_shear    = w_k    .* k_grid.^(-2) .* Pk_grid
weight_gal_shear = w_k   .* Pk_grid
```

Shape: `(96,)` each.

The different \( k \)-power factors \( k^2 \), \( k^{-2} \), and \( k^0 \) arise from the spin-0 (galaxy) and spin-2 (shear) nature of the respective fields, and reflect the prefactors in the non-Limber projection formulae. Specifically:
- **Galaxy–Galaxy** (\( \beta = 2 \)): \( k^2 P(k) \) — the extra \( k^2 \) comes from the fact that the galaxy kernel projects onto scalar overdensity modes.
- **Shear–Shear** (\( \beta = -2 \)): \( k^{-2} P(k) \) — the two factors of \( (\ell+2)!/(\ell-2)! \) in the shear angular power spectrum are compensated by two inverse \( k \)-factors from the relation between the convergence and the lensing potential.
- **Galaxy–Shear** (\( \beta = 0 \)): \( P(k) \) with no \( k \) modification, as the mixed cross-correlation sits between the two pure cases.

---

## 11. Multipole Prefactors \( f_\ell \)

The angular power spectra for spin-0 and spin-2 fields differ by \( \ell \)-dependent geometric factors. The code computes:

```julia
factorial_frac(ℓ) = (ℓ+2)(ℓ+1)ℓ(ℓ-1)
get_ell_prefactor(Galaxy, Galaxy, ℓ)  = (2/π)
get_ell_prefactor(Galaxy, Shear, ℓ)  = (2/π) * sqrt(factorial_frac(ℓ))
get_ell_prefactor(Shear, Shear, ℓ)   = (2/π) * factorial_frac(ℓ)
```

- `pref_gg`: **Galaxy–galaxy prefactor** \( \frac{2}{\pi} \), arising from the standard relation between the angular and 3D power spectra via the spherical harmonic decomposition.
- `pref_ss`: **Shear–shear prefactor** \( \frac{2}{\pi} \cdot (\ell+2)(\ell+1)\ell(\ell-1) \). The combinatorial factor converts spin-2 shear correlations (related to second derivatives of the lensing potential) to convergence power spectra via the relation \( \tilde{\gamma}_\ell = \sqrt{\frac{(\ell+2)!}{(\ell-2)!}} \kappa_\ell \), which gives the square \( (\ell+2)(\ell+1)\ell(\ell-1) \) in the power spectrum.
- `pref_gs`: **Galaxy–shear prefactor** \( \frac{2}{\pi} \cdot \sqrt{(\ell+2)(\ell+1)\ell(\ell-1)} \), the geometric mean of the two.

Stored as 1D arrays of length \( N_\ell = 100 \): `pref_gg`, `pref_ss`, `pref_gs`.

---

## 12. Intermediate Power Spectrum Matrix \( S_\ell(k', k'') \)

Computed via `@tullio` contractions:

```julia
@tullio S_lkk_gg[kp, kpp, li] = pref_gg[li] * weight_gal[k] * W_final_gal[li, k, kp] * W_final_gal[li, k, kpp]
@tullio S_lkk_ss[kp, kpp, li] = pref_ss[li] * weight_shear[k] * W_final_shear[li, k, kp] * W_final_shear[li, k, kpp]
@tullio S_lkk_gs[kp, kpp, li] = pref_gs[li] * weight_gal_shear[k] * W_final_shear[li, k, kp] * W_final_gal[li, k, kpp]
```

Shape: `(96, 96, 100)` = `(N_\chi, N_\chi, N_\ell)`.

The index `k` (summed over, dim = 96) represents the **first \( k \)-integral** weighted by the power spectrum. Indices `kp` and `kpp` label the two remaining comoving-distance dimensions. Physically, \( S_\ell(k', k'') \) is an intermediate quantity that factorizes the three-dimensional radial integration into steps: the first \( k \)-integration (with power spectrum weight) has already been performed, while the remaining two integrations over \( k' \) and \( k'' \) couple to the window function projections. This structure follows from writing the 3D matter power spectrum non-Limber projection as a triple \( k \)-integral:
\[
C_\ell \propto \int dk \int dk' \int dk'' \, P(k) \, W(k') \, W(k'') \, \tilde{W}(k, k', k'')
\]

---

## 13. Final Angular Power Spectra \( S_\ell \)

Computed via:

```julia
weight_kp = w_kp .* kp_grid
@tullio S_l_gg[li] = weight_kp[kp] * S_lkk_gg[kp, kp, li]
@tullio S_l_ss[li] = weight_kp[kp] * weight_kp[kpp] * S_lkk_ss[kp, kpp, li]
@tullio S_l_gs[li] = weight_kp[kp] * weight_kp[kpp] * S_lkk_gs[kp, kpp, li]
```

Shape: `(100,)` = `(N_\ell,)` each.

- `weight_kp = w_kp * kp_grid`: Quadrature weight for the remaining \( k' \) integration, including a factor of \( k' \) from the measure.
- `S_l_gg`: **Galaxy–galaxy angular power spectrum** \( C_\ell^{gg} \). Note that `S_lkk_gg` is contracted only on the diagonal `kp == kpp`, indicating that the galaxy–galaxy projection reduces to a single remaining \( k \)-integral (as opposed to a double integral for the other probes). This is consistent with the structure where the galaxy kernel has an additional \( k \)-power that allows dimensional reduction.
- `S_l_ss`: **Shear–shear angular power spectrum** \( C_\ell^{\kappa\kappa} \), obtained by double contraction over both `kp` and `kpp`.
- `S_l_gs`: **Galaxy–shear cross angular power spectrum** \( C_\ell^{g\kappa} \), similarly double-contracted.

The plotted quantity is \( \ell(\ell+1) C_\ell \), which for a scale-invariant Harrison–Zel'dovich spectrum would be flat, and is the conventional representation that makes the acoustic features visually prominent.

---

## 14. Multipole Array \( \ell \)

- `ℓ = LinRange(2, 200, 100)`: The **angular multipole** \( \ell \), a dimensionless integer (treated as continuous here) that labels the angular Fourier modes on the sphere, analogous to wavenumber \( k \) in 3D. Low \( \ell \) correspond to large-angle modes; \( \ell \sim 100 \) corresponds to angular scales of order a few degrees. The range 2–200 captures the largest scales accessible to wide-field surveys while remaining in the quasi-linear regime where the power spectrum can be reliably predicted.

---

## 15. Summary Table of Key Arrays

| Symbol | Variable name | Shape | Units | Physical meaning |
|---|---|---|---|---|
| \( z \) | `z_range`, `zed` | `(96,)` | dimensionless | Cosmological redshift |
| \( \chi \) | `x_range`, `chi` | `(96,)` | Mpc/h | Comoving distance |
| \( E(z) \) | `Blast.compute_adimensional_hubble_factor` | scalar | dimensionless | Dimensionless Hubble rate |
| \( H(z) \) | `Blast.compute_hubble_factor` | scalar | km/s/Mpc | Physical Hubble rate |
| \( c/H(z) \) | `inv_Hubble_array` | `(96,)` | Mpc/h | Radial comoving line element \( d\chi/dz \) |
| \( D(z) \) | `D_growth_array`, `growth` | `(96,)` | dimensionless | Linear growth factor, normalized to 1 at \( z=0 \) |
| \( n(z) \) | `nz_norm` | `(96,)` | 1 | Normalized source redshift distribution |
| \( b(z) \) | `bz_array`, `bias` | `(96,)` | dimensionless | Linear galaxy bias |
| \( q(z) \) | `lens_int_array` | `(96,)` | 1 | Lensing efficiency integral |
| \( W^g(z) \) | `gal_prefact_W` | `(96,)` | (Mpc/h)² | Galaxy kernel (unnormalized) |
| \( W^\kappa(z) \) | `shear_prefac_W` | `(96,)` | (Mpc/h)⁻¹ | Shear kernel |
| \( c_n^g \) | `cheb_coeff_gal` | `(119,)` | — | Chebyshev coefficients of galaxy kernel |
| \( c_n^\kappa \) | `cheb_coeff_shear` | `(119,)` | — | Chebyshev coefficients of shear kernel |
| \( T_n(k_m) \) | `T` | `(120, N)` | dimensionless | Chebyshev polynomial matrix |
| \( j_\ell(\chi_i k_m) \) | `Bessel` | `(96, N)` | dimensionless | Spherical Bessel function matrix |
| \( \tilde{W}^{(\ell)}_{i,p,n} \) | `W_tilde` | `(100, 96, 96, 119)` | — | Chebyshev–Bessel overlap tensor |
| \( W^{(\ell)}_{\text{final},i,p} \) | `W_final_gal/shear` | `(100, 96, 96)` | — | \( k \)-integrated Bessel transform of kernel |
| \( P(k,z) \) | `Pklin`, `Pknonlin` | `(Nk, Nz)` | (Mpc/h)³ | Linear/nonlinear matter power spectrum |
| \( k \) | `k`, `k_grid` | `(Nk,)` | h/Mpc | Fourier wavenumber |
| \( w_k \) | `w_k` | `(96,)` | — | CC quadrature weights in \( k \)-space |
| \( S_\ell(k',k'') \) | `S_lkk_gg/ss/gs` | `(96, 96, 100)` | — | Intermediate \( k \)-integrated power spectrum |
| \( C_\ell^{gg} \) | `S_l_gg` | `(100,)` | — | Galaxy–galaxy angular power spectrum |
| \( C_\ell^{\kappa\kappa} \) | `S_l_ss` | `(100,)` | — | Shear–shear angular power spectrum |
| \( C_\ell^{g\kappa} \) | `S_l_gs` | `(100,)` | — | Galaxy–shear cross angular power spectrum |
| \( \ell \) | `ℓ` | `(100,)` | dimensionless | Angular multipole |
| \( \beta \) | `β` | scalar | dimensionless | \( k \)-power exponent per probe type |
