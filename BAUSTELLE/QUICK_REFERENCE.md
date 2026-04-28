# Quick Reference: Equations → Code Mapping & Checklists

## Master Equation-to-Code Index

| # | Equation | Mathematical Form | Julia Function | File | Status |
|---|----------|-------------------|-----------------|------|--------|
| **LAYER 0: Cosmological Parameters** |
| 0.1 | Ω_k | $\Omega_k = 1 - \Omega_m - \Omega_\Lambda$ | `CosmoParams()` | cosmo.jl | ✓ Exists |
| 0.2 | H₀ | $H_0 = 100h \text{ km/s/Mpc}$ | `H0(cosmo)` | cosmo.jl | ✓ Exists |
| **LAYER 1: Background Evolution** |
| 1.1 | E(z) | $E(z) = \sqrt{\Omega_m(1+z)^3 + \Omega_k(1+z)^2 + \Omega_\Lambda w(a)}$ | `E(z, cosmo)` | background.jl | ✓ Exists |
| 1.2 | w(a) | $w(a) = w_0 + w_a(1-a)$ | (internal to E) | background.jl | ✓ Exists |
| 1.3 | χ(z) | $\chi(z) = (c/H_0) \int_0^z dz'/E(z')$ | `comoving_distance(z, cosmo)` | background.jl | ⚠ Needs expansion |
| 1.4 | D(z) | ODE from Friedmann equations | `growth_function(z, cosmo)` | background.jl | ⚠ Needs expansion |
| 1.5 | H(z) | $H(z) = H_0 \cdot E(z)$ | (derived from E) | background.jl | ✓ Trivial |
| 1.6 | Background Cache | All of above on grid | `Background` struct | background.jl | ⚠ Needs creation |
| **LAYER 2A: Matter Power Spectrum** |
| 2a.1 | P_lin(k) | Linear power from CLASS/CAMB | `load_linear_power_spectrum()` | projected_matter.jl | ⚠ Needs implementation |
| 2a.2 | P_nl(k) | $P_nl(k) = P_lin(k) \cdot [1 + ...]$ | `nonlinear_power_spectrum()` | projected_matter.jl | ⚠ Needs implementation |
| 2a.3 | P^3D(k,χ,χ') | $P = D(\chi)D(\chi')P_{lin}(k)$ | `power_spectrum_3d()` | projected_matter.jl | ⚠ Needs implementation |
| **LAYER 2B: Spherical Bessel** |
| 2b.1 | j_ℓ(x) | $j_\ell(x) = \sqrt{\pi/2x} J_{\ell+1/2}(x)$ | `spherical_bessel(ell, x)` | spherical_bessel.jl | ❌ NEW FILE |
| 2b.2 | j'_ℓ(x) | $dj_\ell/dx = j_{\ell-1} - (\ell+1)j_\ell/x$ | `spherical_bessel_derivative()` | spherical_bessel.jl | ❌ NEW FILE |
| 2b.3 | x_{ℓ,n} | Zeros where $j_\ell(x) = 0$ | `spherical_bessel_zeros()` | spherical_bessel.jl | ❌ NEW FILE |
| 2b.4 | Zeros Cache | Precomputed zeros | `get_spherical_bessel_zeros()` | spherical_bessel.jl | ❌ NEW FILE |
| **LAYER 3: Window Functions** |
| 3.1 | f^den(χ) | $f = (H/c) b(z) n(z) D(z)$ | `density_weight_function()` | sfb_kernels.jl | ❌ NEW FILE |
| 3.2 | f^lens(χ) | $f = (3/2)\Omega_m(H_0/c)^2 \int ...$ | `lensing_weight_function()` | sfb_kernels.jl | ❌ NEW FILE |
| 3.3 | W_ℓ^A(k,χ) | $W = f^A(\chi) j_\ell(k\chi)$ | `window_function()` | sfb_kernels.jl | ❌ NEW FILE |
| 3.4 | W̃_ℓ^A(k,k') | $\widetilde{W} = \int d\chi \chi^2 W j_\ell(k'\chi)$ | `hankel_transform()` | sfb_kernels.jl | ❌ NEW FILE |
| 3.5 | W̃_ℓ grid | Hankel over k' array | `hankel_transform_grid()` | sfb_kernels.jl | ❌ NEW FILE |
| **LAYER 4: 3D Correlations** |
| 4.1 | S_ℓ^AB(k₁,k₂) | $S = N_\ell \int dk k^2 P W̃^A W̃^B$ | `compute_correlation_function()` | sfb_decomposition.jl | ❌ NEW FILE |
| 4.2 | S_ℓ grid | Correlation on (k₁,k₂) mesh | `compute_correlation_grid()` | sfb_decomposition.jl | ❌ NEW FILE |
| 4.3 | SFB Cache | Precomputed Hankel transforms | `SFBCorrelationCache` struct | sfb_decomposition.jl | ❌ NEW FILE |
| **LAYER 5: Final Observables** |
| 5.1 | C_ℓ^AB(z_i,z_j) | Integrate S_ℓ over wavenumbers | `angular_power_spectrum()` | integrals.jl | ⚠ Needs expansion |
| 5.2 | Validation | Compare vs Limber | `validate_against_limber()` | integrals.jl | ❌ NEW FUNCTION |

---

## File-by-File Implementation Checklist

### **cosmo.jl** (Existing - Minor Extensions)

**Status**: Mostly complete, may need small additions

```
[ ] Struct CosmoParams defined
    [ ] Ω_m, Ω_Λ, Ω_k
    [ ] h, w0, wa
    [ ] σ8, n_s (if needed)
    [ ] Validation in constructor
    
[ ] Derived quantities:
    [ ] H0(cosmo) → H₀ in km/s/Mpc
    [ ] H0_Mpc(cosmo) → H₀ in 1/Mpc
    [ ] validate_consistency() → check Ω_k = 1 - Ω_m - Ω_Λ
```

---

### **background.jl** (Existing - Major Expansion)

**Status**: Some functions exist, need systematic completion

```
✓ EXISTING:
[ ] E(z, cosmo) function exists

⚠ NEEDS IMPLEMENTATION:
[ ] comoving_distance(z, cosmo)
    [ ] Implement numerical integration (quadgk)
    [ ] Handle z ≈ 0 case
    [ ] Test against precomputed data/background/chi.npy
    
[ ] growth_function(z, cosmo)
    [ ] Set up 2nd-order ODE
    [ ] Integrate z → 0
    [ ] Normalize D(z=0) = 1
    [ ] Test against reference values
    
[ ] Hubble parameter H(z, cosmo)
    [ ] H(z) = H₀ * E(z)
    [ ] Return in 1/Mpc (not km/s/Mpc)

[ ] Background struct
    [ ] Store z, E(z), χ(z), D(z), H(z) on grid
    [ ] Create interpolators for fast lookup
    [ ] Include inverse interpolators (z ↔ χ)

[ ] Precomputation function
    [ ] precompute_background(z_max, n_pts, cosmo)
    [ ] Cache with hash of cosmo params
    [ ] Save to data/ if needed

TESTS:
[ ] @test E(0) ≈ 1
[ ] @test χ(z) monotonically increasing
[ ] @test D(z=0) ≈ 1
[ ] @test H(z) > 0 for all z
[ ] Compare vs data/background/chi.npy (error < 0.1%)
```

---

### **projected_matter.jl** (Existing - Major Expansion)

**Status**: Minimal, needs complete implementation

```
❌ NEEDS IMPLEMENTATION:

[ ] load_linear_power_spectrum(filename)
    [ ] Parse CLASS/CAMB output
    [ ] Create log-linear interpolator
    [ ] Test interpolation accuracy
    [ ] Handle boundaries (extrapolation)
    
[ ] nonlinear_power_spectrum(k, z, P_lin, cosmo)
    [ ] Implement HaloFit or similar
    [ ] Scale P_lin by D²(z)
    [ ] Test against reference codes
    
[ ] power_spectrum_3d(k, χ, χ', P_lin, background, cosmo)
    [ ] Convert χ → z
    [ ] Get D(z) from background
    [ ] Multiply: P = D(z) * D(z') * P_lin(k)
    [ ] Vectorize over arrays
    
[ ] Integration setup
    [ ] Define k grid for integrals
    [ ] Precompute P_lin on grid
    [ ] Setup caching

TESTS:
[ ] P_lin(k) > 0 for all k
[ ] P_lin monotonically decreases at high k
[ ] P_nl ≤ P_lin (approximately)
[ ] Comparison with HaloFit reference code
```

---

### **spherical_bessel.jl** (NEW FILE)

**Status**: Must create from scratch

```
❌ CREATE NEW FILE:

[ ] spherical_bessel(ell, x)
    [ ] Import SpecialFunctions.jl
    [ ] Handle x ≈ 0 case: j_0(0)=1, j_ℓ(0)=0
    [ ] Use besselj(ell+0.5, x) * sqrt(π/(2x))
    [ ] Vectorized version with @.
    [ ] Benchmark vs naive implementation
    
[ ] spherical_bessel_derivative(ell, x)
    [ ] Implement recurrence: j'_ℓ = j_{ℓ-1} - (ℓ+1)j_ℓ/x
    [ ] Vector version
    [ ] Test against finite differences
    
[ ] spherical_bessel_zeros(ell, n_zeros; x_max=1000)
    [ ] Find zeros using Roots.jl
    [ ] Bracket zeros using Bessel oscillation properties
    [ ] Return sorted array
    [ ] Handle errors gracefully
    
[ ] get_spherical_bessel_zeros(ell, n_zeros)
    [ ] Cache results in BESSEL_ZEROS_CACHE
    [ ] Retrieve from cache on repeat calls
    [ ] Thread-safe (if needed)

VALIDATION:
[ ] j_0(π) ≈ 0 (first zero)
[ ] j_1(π) ≠ 0 (not a zero of j_0)
[ ] Zeros match GSL reference values
[ ] Derivative matches finite-difference test
[ ] Benchmark: <1 ms per Bessel call
```

---

### **sfb_kernels.jl** (NEW FILE)

**Status**: Must create from scratch

```
❌ CREATE NEW FILE:

[ ] SurveyParams struct
    [ ] Store n(z), b(z) as functions
    [ ] Loader from file
    [ ] Getters: get_number_density(), get_bias()
    
[ ] density_weight_function(chi, background, survey)
    [ ] f^den(χ) = (H/c) * b(z) * n(z) * D(z)
    [ ] Convert χ → z
    [ ] Retrieve all components from background
    [ ] Handle numerical stability
    
[ ] lensing_weight_function(chi, background, survey)
    [ ] f^lens(χ) = (3/2)*Ω_m*(H0/c)² * ∫ dχ' n(χ') W^lens
    [ ] More complex: involves 2nd integral
    [ ] (Optional: implement if needed)
    
[ ] window_function(ell, k, chi, f_weight)
    [ ] W_ℓ(k,χ) = f(χ) * j_ℓ(k*χ)
    [ ] Calls spherical_bessel()
    [ ] Simple, fast function
    
[ ] hankel_transform(ell, k, k1, z_min, z_max, f_weight, ...)
    [ ] Define integrand: χ² * f(χ) * j_ℓ(kχ) * j_ℓ(k1χ)
    [ ] Use quadgk for adaptive integration
    [ ] Set rtol=1e-5
    [ ] Return scalar
    
[ ] hankel_transform_grid(ell, k, k1_arr, ...)
    [ ] Vectorized over k1_arr
    [ ] More efficient than loop of hankel_transform
    [ ] Return vector [n_k1]

VALIDATION:
[ ] W_ℓ = 0 where f = 0
[ ] Hankel(const weight) matches analytical form
[ ] Convergence: error decreases with integration tolerance
[ ] Benchmark: <100 ms per Hankel transform
[ ] Compare vs Limber approximation for large ℓ
```

---

### **sfb_decomposition.jl** (NEW FILE)

**Status**: Must create from scratch

```
❌ CREATE NEW FILE:

[ ] SFBCorrelationCache struct
    [ ] ell::Int
    [ ] k_vals, k1_vals, k2_vals arrays
    [ ] W_tilde_A, W_tilde_B matrices
    [ ] P_lin_vals vector
    [ ] Interpolators (W_A_itp, W_B_itp, P_itp)
    
[ ] SFBCorrelationCache(ell, k_min, k_max, n_k, ...)
    [ ] Constructor that precomputes everything
    [ ] Setup logarithmic k-grid
    [ ] Compute Hankel transforms (expensive!)
    [ ] Create interpolators
    [ ] Store in cache
    
[ ] compute_correlation_function(ell, k1, k2, cache)
    [ ] Define integrand: k² * P_lin(k) * W̃^A(k,k1) * W̃^B(k,k2)
    [ ] Interpolate W̃ at (k, k1) and (k, k2)
    [ ] Use quadgk over k
    [ ] Return S_ℓ(k1, k2)
    
[ ] compute_correlation_grid(ell, k1_arr, k2_arr, cache)
    [ ] Vectorized: nested loop over k1, k2
    [ ] Return matrix [n_k1 × n_k2]
    [ ] Parallelize if needed
    
[ ] Precomputation wrapper
    [ ] precompute_sfb_cache(ell_arr, survey, background, P_lin, cosmo)
    [ ] Return Dict[ell] → SFBCorrelationCache
    [ ] Handle all ℓ values needed

VALIDATION:
[ ] S_ℓ(k1, k2) symmetric in k1, k2 (if W̃^A = W̃^B)
[ ] S_ℓ > 0 everywhere
[ ] Compare with Limber reference (error < 1% for large ℓ)
[ ] Convergence: error vs number of k-grid points
[ ] Benchmark: <1 s for full grid computation
```

---

### **integrals.jl** (Existing - Major Expansion)

**Status**: Some Limber code exists, need SFB extension

```
✓ EXISTING:
[ ] Limber-based computations (reference)

⚠ NEEDS IMPLEMENTATION:

[ ] angular_power_spectrum(ell, z_i, z_j, sfb_cache, cosmo)
    [ ] Integrate S_ℓ over three wavenumbers
    [ ] Use hcubature (3D quadrature) or similar
    [ ] Set appropriate tolerances
    [ ] Return C_ℓ^AB(z_i, z_j)
    
[ ] angular_power_spectrum_grid(ell_arr, z_bins, sfb_cache, cosmo)
    [ ] Vectorized: compute C_ℓ for all ell and redshift pairs
    [ ] Return tensor [n_ell × n_z_bins × n_z_bins]
    
[ ] validate_against_limber(ell_arr, z_bins, C_sfb, C_limber)
    [ ] Load precomputed Limber data
    [ ] Compute relative errors
    [ ] Generate comparison plots
    [ ] Document convergence behavior
    
[ ] Convergence analysis
    [ ] Error vs number of SFB modes
    [ ] Error vs integration tolerance
    [ ] Error vs k-grid resolution

TESTS:
[ ] C_ℓ^AB symmetric for A=B
[ ] C_ℓ > 0 (acoustic oscillations can be negative, but overall positive)
[ ] Decreases with ℓ (roughly)
[ ] Comparison with Limber: error < 1%
[ ] Scaling: time ∝ n_ell * n_zbins²
```

---

## Implementation Timeline & Milestones

### **Week 1-2: Layers 0-1 (Background)**
```
MON: Review existing cosmo.jl and background.jl code
TUE: Implement comoving_distance(z) + test
WED: Implement growth_function(z) + test
THU: Create Background struct + precompute_background()
FRI: Validation notebook: plot E(z), χ(z), D(z), compare with data/background/
```

### **Week 2-3: Layer 2A (Power Spectrum)**
```
MON: load_linear_power_spectrum() + test
TUE: nonlinear_power_spectrum() + HaloFit implementation
WED: power_spectrum_3d() + vectorization
THU: Benchmark and optimization
FRI: Validation notebook: plot P_lin, P_nl vs reference
```

### **Week 3-4: Layer 2B (Spherical Bessel)**
```
MON: Create spherical_bessel.jl file
TUE: spherical_bessel(ell, x) + vectors
WED: spherical_bessel_zeros() with caching
THU: Benchmark vs GSL/scipy, optimize
FRI: Validation notebook: Bessel functions, zeros, derivatives
```

### **Week 4-6: Layer 3 (Window Functions)**
```
MON: Create sfb_kernels.jl file
TUE: density_weight_function() + SurveyParams
WED: window_function() + test
THU: hankel_transform() single call
FRI: hankel_transform_grid() vectorized

MON: Benchmark + optimize quadrature
TUE: Cache strategies
WED: Validation notebook: compare Hankel transforms
THU: Integration with Layer 2B (Bessel)
FRI: Performance profiling
```

### **Week 6-7: Layer 4 (Correlations)**
```
MON: Create sfb_decomposition.jl
TUE: SFBCorrelationCache struct + constructor
WED: compute_correlation_function(ell, k1, k2)
THU: compute_correlation_grid() + vectorization
FRI: Validation notebook: S_ℓ(k1, k2) heatmaps, vs Limber
```

### **Week 7-8: Layer 5 (Final Observables)**
```
MON: angular_power_spectrum() implementation
TUE: angular_power_spectrum_grid() for full output
WED: validate_against_limber() + convergence analysis
THU: Comprehensive validation notebook
FRI: Documentation + code review
```

### **Week 8+: Optimization & Production**
```
Use @turbo, @tullio for performance
Parallelize where possible
Add error propagation
Create publication-ready code
```

---

## Quick Lookup: "Where Does This Go?"

**Q: I want to implement the growth factor D(z) ODE**  
A: → `src/background.jl`, function `growth_function(z::Float64, cosmo::CosmoParams)`

**Q: I need to compute spherical Bessel function j_2(3.5)**  
A: → `src/spherical_bessel.jl`, function `spherical_bessel(2, 3.5)`

**Q: Where do I handle the survey number density n(z)?**  
A: → `src/sfb_kernels.jl`, struct `SurveyParams`, used in `density_weight_function()`

**Q: How do I compute the Hankel transform W̃_ℓ(k, k')?**  
A: → `src/sfb_kernels.jl`, function `hankel_transform(ell, k, k_prime, ...)`

**Q: Where is the final angular power spectrum C_ℓ computed?**  
A: → `src/integrals.jl`, function `angular_power_spectrum(ell, z_i, z_j, ...)`

**Q: Which file has all the precomputation caching?**  
A: → Multiple files:
  - `background.jl`: Background struct
  - `spherical_bessel.jl`: Bessel zeros cache
  - `sfb_decomposition.jl`: SFBCorrelationCache

**Q: Where do I validate against Limber reference?**  
A: → `src/integrals.jl`, function `validate_against_limber()`, + notebook

---

## Running Full Pipeline

Once all pieces are implemented:

```julia
# 1. Setup
cosmo = CosmoParams(Ω_m=0.3, Ω_Λ=0.7, h=0.67)
background = precompute_background(10.0, 500, cosmo)

# 2. Load data
P_lin = load_linear_power_spectrum("data/pk.npy")
survey = SurveyParams(n_z_file, b_z_file)

# 3. Precompute SFB
ell_arr = [0, 2, 4, 10, 50, 100]
sfb_cache = precompute_sfb_cache(ell_arr, survey, background, P_lin, cosmo)

# 4. Compute observables
z_bins = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.5)]
C_ell = angular_power_spectrum_grid(ell_arr, z_bins, sfb_cache, cosmo)

# 5. Validate
C_limber = load_limber_reference()
residuals = validate_against_limber(ell_arr, z_bins, C_ell, C_limber)
println("Max relative error: $(residuals["max_rel_error"])")
```

---

## Debugging Tips

| Problem | Solution | File |
|---------|----------|------|
| `chi(z)` diverges | Check E(z) → quadgk convergence | background.jl |
| `D(z)` has kink | ODE solver stiffness, use RadauIIA | background.jl |
| `j_ℓ(x)` wrong values | Check Bessel vs reference (scipy) | spherical_bessel.jl |
| `hankel_transform` very slow | Reduce rtol, or use adaptive grid | sfb_kernels.jl |
| `S_ℓ(k1, k2)` has NaN | Check W̃ interpolation out-of-range | sfb_decomposition.jl |
| `C_ℓ` doesn't match Limber | Compare ℓ-dependence, check normalization | integrals.jl |

---

## Key Paper References

- **BLAST paper**: Chiarenza et al. 2024 (arXiv:2410.03632)
  - Equations 15-20: Core formalism
  - Figure 3: Accuracy vs Limber
  
- **Equations.tex** (your file):
  - Generalized BLAST form (Eq. 1)
  - Hankel transform definition (Eq. 2-3)
  - Density window function (Eq. 4)

- **Related**: Limber approximation limits, growth rate formalism, weak lensing kernels

