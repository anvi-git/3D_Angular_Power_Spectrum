# Complete Implementation Guide: Spherical Fourier-Bessel Decomposition for 3D Power Spectra

**Project**: Beyond Limber Angular Power Spectra via Spherical Fourier-Bessel (SFB) Decomposition  
**Framework**: BLAST (arXiv:2410.03632v1)  
**Language**: Julia  
**Timeline**: Incremental development with validation at each layer

---

## Table of Contents

1. [Overview: What We're Building](#overview)
2. [Dependency Chain: What Depends on What](#dependency-chain)
3. [Layer-by-Layer Implementation](#layer-by-layer-implementation)
4. [Mathematical Foundations](#mathematical-foundations)
5. [Code Structure & File Mapping](#code-structure--file-mapping)
6. [Detailed Implementation Steps](#detailed-implementation-steps)

---

## Overview

### The Problem We're Solving

From `equations.tex`, we start with the **generalized BLAST formula** for observable correlations:

$$C_\ell^{AB}(z_i, z_j) = N_\ell^{AB} \int d\chi \int d\chi' \int dk \, k^2 \, P^{AB}(k,\chi,\chi') \, j_\ell(k\chi) \, j_\ell(k\chi')$$

**The challenge**: This triple integral is computationally expensive. BLAST solves it by:

1. **Using Chebyshev polynomials** for efficient precomputation (current BLAST method)
2. **Using Spherical Fourier-Bessel (SFB) decomposition** (your extension) to achieve full 3D treatment beyond Limber

### Key Innovation: Hankel Transform Strategy

Instead of computing the triple integral directly, we decompose using **Hankel transforms**:

$$\widetilde{W}_\ell^A(k, k_1) = \int_0^\infty d\chi \, \chi^2 \, W_\ell^A(k,\chi) \, j_\ell(k_1 \chi)$$

This transforms the 3D integral into **separable 1D integrals** over redshift space, which is dramatically more efficient.

---

## Dependency Chain: What Depends on What

```
┌─────────────────────────────────────────────────────────────────┐
│  LAYER 0: COSMOLOGICAL PARAMETERS (cosmo.jl)                   │
│  ├─ Cosmological constants (Ω_m, Ω_Λ, h, w0, wa, ...)          │
│  └─ Redshift range and number of bins                           │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│  LAYER 1: BACKGROUND EVOLUTION (background.jl)                 │
│  ├─ E(z) = H(z)/H₀ [depends on cosmo params]                   │
│  ├─ χ(z) = comoving distance [depends on E(z)]                 │
│  ├─ D(z) = growth function [depends on cosmo params]           │
│  └─ H(z) = Hubble parameter [depends on E(z)]                  │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│  LAYER 2A: MATTER POWER SPECTRUM (projected_matter.jl)         │
│  ├─ P_lin(k) = linear power spectrum [from CLASS or CAMB]      │
│  ├─ P_nl(k) = nonlinear power spectrum [from HaloFit]         │
│  └─ P_3D(k, χ, χ') = growth weighted power [depends on D(z)]  │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│  LAYER 2B: SPHERICAL BESSEL FUNCTIONS (spherical_bessel.jl)    │
│  ├─ j_ℓ(x) = spherical Bessel functions                        │
│  ├─ j'_ℓ(x) = derivatives                                      │
│  ├─ zeros(ℓ, n_zeros) = zeros of j_ℓ [used for SFB grid]      │
│  └─ recurrence relations for efficiency                        │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│  LAYER 3: WINDOW FUNCTIONS & HANKEL TRANSFORMS (sfb_kernels.jl)│
│  ├─ f^A(χ) = weight function [depends on Background + probe]   │
│  │   - f^den(χ) = H(χ)·b(χ)·n(χ)·D(χ)/c                       │
│  │   - f^lens(χ) = ∫ dχ' n(χ') W^lens(χ,χ')/χ                │
│  ├─ W^A_ℓ(k,χ) = f^A(χ)·j_ℓ(k·r(χ))                          │
│  └─ W̃^A_ℓ(k,k₁) = Hankel transform of W^A [integrates χ]     │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│  LAYER 4: SFB COEFFICIENTS (sfb_decomposition.jl)              │
│  ├─ Compute W̃ on full redshift grid [calls sfb_kernels]       │
│  ├─ SFB basis expansion a_νℓm = projection of P onto basis    │
│  └─ Store coefficients for reuse                               │
└─────────────────────┬───────────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────────┐
│  LAYER 5: FINAL OBSERVABLES (integrals.jl, sfb_decomposition.jl)
│  ├─ S^AB_ℓ(k₁, k₂) = spherical correlation function           │
│  ├─ C^AB_ℓ = integrated 3x2pt angular power spectrum           │
│  └─ Validation: Compare against Limber reference               │
└─────────────────────────────────────────────────────────────────┘
```

### Conceptual Flow

```
Cosmological Model (z, Ω_m, ...)
         ↓
Background Evolution (E(z), χ(z), D(z), ...)
         ↓
Power Spectrum P(k, z)
         ↓
Window Functions f^A(χ), f^B(χ)
         ↓
Spherical Bessel Functions j_ℓ(x)
         ↓
Hankel Transforms: W̃_ℓ(k, k')
         ↓
3D Correlations S_ℓ(k₁, k₂)
         ↓
Angular Power Spectrum C_ℓ (final product)
```

---

## Layer-by-Layer Implementation

### LAYER 0: Cosmological Parameters (`cosmo.jl`)

**Purpose**: Define and store the cosmological model.

**Mathematical quantities**:
- Ω_m, Ω_Λ, Ω_k: matter, dark energy, curvature density fractions
- h = H₀ / (100 km/s/Mpc): normalized Hubble constant
- w₀, wₐ: dark energy equation of state parameters
- σ₈: matter clustering amplitude
- n_s: spectral index

**Pseudocode**:
```julia
struct CosmoParams
    Ω_m::Float64      # Matter density
    Ω_Λ::Float64      # Dark energy density
    Ω_k::Float64      # Curvature (= 1 - Ω_m - Ω_Λ)
    h::Float64        # H₀ / 100
    w0::Float64       # w₀ (dark energy EoS today)
    wa::Float64       # wₐ (dark energy EoS evolution)
    σ8::Float64       # Cluster amplitude
    ns::Float64       # Spectral index
end

function CosmoParams(; Ω_m=0.3, Ω_Λ=0.7, h=0.67, w0=-1.0, wa=0.0, σ8=0.8, ns=0.96)
    # Validate: Ω_m + Ω_Λ + Ω_k = 1
    Ω_k = 1.0 - Ω_m - Ω_Λ
    return CosmoParams(Ω_m, Ω_Λ, Ω_k, h, w0, wa, σ8, ns)
end

# Derived quantities
H₀(cosmo::CosmoParams) = 100 * cosmo.h  # km/s/Mpc
H₀_Mpc(cosmo::CosmoParams) = H₀(cosmo) / 299792.458  # 1/Mpc
```

**File location**: `src/cosmo.jl` (already exists, extend if needed)

---

### LAYER 1: Background Evolution (`background.jl`)

**Purpose**: Compute cosmological background quantities that evolve with redshift.

**Key Functions to Implement**:

#### 1a. **Hubble Parameter E(z)**

From Friedmann equation:
$$E(z) = \frac{H(z)}{H_0} = \sqrt{\Omega_m (1+z)^3 + \Omega_k(1+z)^2 + \Omega_\Lambda \, w(z)}$$

where for w₀-wₐ parameterization:
$$w(a) = w_0 + w_a(1-a) = w_0 + w_a \frac{z}{1+z}$$

**Implementation**:
```julia
"""
    E(z::Float64, cosmo::CosmoParams)::Float64

Hubble parameter normalized to H₀: H(z) = H₀ * E(z)
"""
function E(z::Float64, cosmo::CosmoParams)::Float64
    a = 1.0 / (1.0 + z)
    w_z = cosmo.w0 + cosmo.wa * (1.0 - a)
    
    term_m = cosmo.Ω_m * (1.0 + z)^3
    term_k = cosmo.Ω_k * (1.0 + z)^2
    term_Λ = cosmo.Ω_Λ * a^(3 * (1 + w_z))
    
    return sqrt(term_m + term_k + term_Λ)
end

# Vectorized version for redshift arrays
E(z_arr::Vector{Float64}, cosmo::CosmoParams) = @. E(z_arr, cosmo)
```

#### 1b. **Comoving Distance χ(z)**

$$\chi(z) = \frac{c}{H_0} \int_0^z \frac{dz'}{E(z')}$$

**Implementation** (using adaptive quadrature):
```julia
"""
    comoving_distance(z::Float64, cosmo::CosmoParams)::Float64

Comoving distance from observer to redshift z (in Mpc)
"""
function comoving_distance(z::Float64, cosmo::CosmoParams)::Float64
    if z < 1e-6
        return 0.0
    end
    
    c_over_H0 = 2997.92458  # Mpc (speed of light / H₀ in Mpc)
    integrand(z_prime) = 1.0 / E(z_prime, cosmo)
    
    integral, _ = quadgk(integrand, 0.0, z, rtol=1e-6)
    return c_over_H0 * integral
end

# For efficiency, precompute on a grid and interpolate
"""
    setup_distance_interpolation(z_max::Float64, n_pts::Int, cosmo::CosmoParams)

Precompute χ(z) on a grid and create fast interpolator
"""
function setup_distance_interpolation(z_max::Float64, n_pts::Int, cosmo::CosmoParams)
    z_grid = range(0, z_max, length=n_pts)
    chi_grid = comoving_distance.(z_grid, Ref(cosmo))
    
    chi_interp = linear_interpolation(z_grid, chi_grid, extrapolation_bc=Line())
    return chi_interp
end
```

#### 1c. **Growth Function D(z)**

$$\frac{d^2D}{da^2} + \frac{3}{a} \frac{dD}{da} + \frac{1}{a^2 E^2(a)} \left(\Omega_m(a) - 2 - \frac{3}{2} \frac{\Omega_\Lambda(a)}{E^2(a)}\right) D = 0$$

With normalized condition: $D(z=0) = 1$.

**Implementation**:
```julia
"""
    growth_function(z::Float64, cosmo::CosmoParams)::Float64

Linear growth function D(z), normalized to D(z=0) = 1
Solves second-order ODE from z to z=0
"""
function growth_function(z::Float64, cosmo::CosmoParams)::Float64
    # For flat ΛCDM with constant w, can use analytic approximation
    # For general case, solve ODE
    
    function growth_ode!(du, u, p, a_vals)
        z_val = 1.0 / a_vals - 1.0
        D, dD_da = u
        
        E_z = E(z_val, cosmo)
        Ω_m_z = cosmo.Ω_m * (1 + z_val)^3 / E_z^2
        
        d2D_da2 = -3.0/a_vals * dD_da - (Ω_m_z / (2.0 * a_vals^2 * E_z^2)) * D
        
        return [dD_da, d2D_da2]
    end
    
    a_init = 1.0e-3  # Start at high z
    a_fin = 1.0 / (1 + z)
    
    prob = ODEProblem(growth_ode!, [1e-3, 0.0], [a_init, a_fin], cosmo)
    sol = solve(prob, RK45())
    
    return sol(a_fin)[1]
end
```

#### 1d. **Precomputation Strategy**

For efficiency, **precompute and cache**:
```julia
"""
    struct Background
        z::Vector{Float64}
        E_z::Vector{Float64}
        chi_z::Vector{Float64}
        D_z::Vector{Float64}
        H_z::Vector{Float64}  # H(z) in 1/Mpc
        interp_E::Interpolation
        interp_chi::Interpolation
        interp_D::Interpolation
        interp_H::Interpolation
    end
"""
```

**File location**: `src/background.jl` (expand existing file)

---

### LAYER 2A: Matter Power Spectrum (`projected_matter.jl`)

**Purpose**: Provide the 3D matter power spectrum $P(k, z)$, which is the core ingredient.

#### 2a-i. **Linear Power Spectrum $P_{lin}(k)$**

**Source**: CLASS or CAMB (precomputed, loaded from file)

**Implementation**:
```julia
"""
    load_linear_power_spectrum(filename::String)::Interpolation

Load precomputed linear power spectrum from file (CLASS/CAMB output)
Returns interpolator P_lin(k) valid for k ∈ [k_min, k_max]
"""
function load_linear_power_spectrum(filename::String)
    data = readdlm(filename)
    k_arr = data[:, 1]      # [1/Mpc]
    P_arr = data[:, 2]      # [Mpc³]
    
    P_interp = log_linear_interpolation(log.(k_arr), log.(P_arr))
    return x -> exp(P_interp(log(x)))
end

# Or: use CLASS.jl directly if installed
# function compute_linear_power_spectrum(cosmo::CosmoParams, k_arr::Vector)::Vector
#     ...
# end
```

#### 2a-ii. **Nonlinear Power Spectrum $P_{nl}(k)$**

**Method**: Halofit (or other phenomenological model)

$$P_{nl}(k) = \frac{P_{lin}(k)}{1 + k \, \sigma_8}^2 \cdot \text{(HaloFit correction)}$$

```julia
"""
    nonlinear_power_spectrum(k::Float64, z::Float64, 
                            P_lin::Function, cosmo::CosmoParams)::Float64

Nonlinear power spectrum using Halofit approximation
"""
function nonlinear_power_spectrum(k::Float64, z::Float64, 
                                  P_lin::Function, cosmo::CosmoParams)::Float64
    D_z = growth_function(z, cosmo)
    D_0 = 1.0  # D(z=0)
    
    # Scale P_lin to redshift z
    P_lin_z = P_lin(k) * (D_z / D_0)^2
    
    # Apply Halofit-like correction (simplified)
    k_nl = 1.0  # Non-linear scale
    correction = 1.0 + (k / k_nl)^2
    
    return P_lin_z / correction^2
end
```

#### 2a-iii. **3D Power Spectrum $P^{AB}(k, \chi, \chi')$**

The full 3D power spectrum includes growth factor weighting:

$$P^{AB}(k, \chi, \chi') = D(\chi) \, D(\chi') \, P_{lin}^{AB}(k)$$

**Implementation**:
```julia
"""
    power_spectrum_3d(k::Float64, chi::Float64, chi_prime::Float64,
                     P_lin::Function, background::Background, 
                     cosmo::CosmoParams, spec_type::String="linear")::Float64

Full 3D power spectrum with growth weighting
spec_type ∈ {"linear", "nonlinear"}
"""
function power_spectrum_3d(k::Float64, chi::Float64, chi_prime::Float64,
                          P_lin::Function, background::Background, 
                          cosmo::CosmoParams, spec_type::String="linear")::Float64
    
    # Convert comoving distance to redshift
    z = background.interp_chi_inv(chi)
    z_prime = background.interp_chi_inv(chi_prime)
    
    # Get growth functions
    D_z = background.interp_D(z)
    D_z_prime = background.interp_D(z_prime)
    
    # Get power spectrum at k
    if spec_type == "linear"
        P_k = P_lin(k)
    else
        P_k = nonlinear_power_spectrum(k, (z + z_prime) / 2, P_lin, cosmo)
    end
    
    return D_z * D_z_prime * P_k
end
```

**File location**: `src/projected_matter.jl` (expand existing file)

---

### LAYER 2B: Spherical Bessel Functions (`spherical_bessel.jl` - NEW FILE)

**Purpose**: Efficiently compute spherical Bessel functions $j_\ell(x)$ and their zeros.

#### 2b-i. **Spherical Bessel Functions $j_\ell(x)$**

$$j_\ell(x) = \sqrt{\frac{\pi}{2x}} J_{\ell+1/2}(x)$$

**Implementation** (using `SpecialFunctions.jl`):
```julia
"""
    spherical_bessel(ell::Int, x::Float64)::Float64

Spherical Bessel function j_ℓ(x)
Uses recurrence relation for efficiency
"""
function spherical_bessel(ell::Int, x::Float64)::Float64
    if abs(x) < 1e-10
        return (ell == 0) ? 1.0 : 0.0
    end
    
    # Use SpecialFunctions.jl
    return besselj(ell + 0.5, x) * sqrt(π / (2x))
end

# Vectorized
spherical_bessel(ell::Int, x_arr::Vector{Float64}) = @. spherical_bessel(ell, x_arr)
```

#### 2b-ii. **Derivative of Spherical Bessel Functions**

$$\frac{dj_\ell}{dx}(x) = j_{\ell-1}(x) - \frac{\ell+1}{x} j_\ell(x)$$

```julia
"""
    spherical_bessel_derivative(ell::Int, x::Float64)::Float64

Derivative of j_ℓ(x) with respect to x
"""
function spherical_bessel_derivative(ell::Int, x::Float64)::Float64
    if abs(x) < 1e-10
        return 0.0
    end
    
    j_ell = spherical_bessel(ell, x)
    j_ell_minus = spherical_bessel(ell - 1, x)
    
    return j_ell_minus - (ell + 1.0) / x * j_ell
end
```

#### 2b-iii. **Zeros of Spherical Bessel Functions**

For SFB decomposition, need zeros $x_{\ell,n}$ where $j_\ell(x_{\ell,n}) = 0$:

```julia
"""
    spherical_bessel_zeros(ell::Int, n_zeros::Int; x_max::Float64=1000.0)::Vector{Float64}

Find first n_zeros zeros of j_ℓ(x) for given ℓ
Returns vector of zeros: [x_ℓ,₁, x_ℓ,₂, ...]
"""
function spherical_bessel_zeros(ell::Int, n_zeros::Int; x_max::Float64=1000.0)::Vector{Float64}
    zeros_ell = Float64[]
    
    # Search for zeros using root-finding
    x_current = ell + 1.0  # Approximate first zero
    
    for n in 1:n_zeros
        # Find zero near x_current using bisection or Newton's method
        f(x) = spherical_bessel(ell, x)
        f_prime(x) = spherical_bessel_derivative(ell, x)
        
        try
            zero_found = find_zero(f, (x_current, x_current + π), 
                                   Bisection())
            push!(zeros_ell, zero_found)
            x_current = zero_found + π * 0.5  # Estimate next zero location
        catch
            @warn "Could not find zero $n for ℓ=$ell"
            break
        end
    end
    
    return zeros_ell
end

# Cache for efficiency
const SPHERICAL_BESSEL_ZEROS_CACHE = Dict()

function get_spherical_bessel_zeros(ell::Int, n_zeros::Int)::Vector{Float64}
    key = (ell, n_zeros)
    if !haskey(SPHERICAL_BESSEL_ZEROS_CACHE, key)
        SPHERICAL_BESSEL_ZEROS_CACHE[key] = spherical_bessel_zeros(ell, n_zeros)
    end
    return SPHERICAL_BESSEL_ZEROS_CACHE[key]
end
```

#### 2b-iv. **Precomputation Strategy**

```julia
"""
    struct SphericalBesselBasis
        ell::Int                     # Multipole order
        n_modes::Int                 # Number of radial modes
        zeros::Vector{Float64}       # Zeros x_ℓ,n
        radii::Vector{Float64}       # r_ℓ,n = x_ℓ,n / k_max (normalized)
    end
"""
```

**File location**: `src/spherical_bessel.jl` (NEW FILE)

---

### LAYER 3: Window Functions & Hankel Transforms (`sfb_kernels.jl` - NEW FILE)

**Purpose**: Compute window functions and their Hankel transforms, which encode the redshift and angular information.

#### 3-i. **Window Function Components**

For **density** probe (galaxies):
$$f^{\text{den}}(\chi) = \frac{H(\chi)}{c} \cdot b(\chi) \cdot n(\chi) \cdot D(\chi)$$

where:
- $H(\chi)$ = Hubble parameter (from comoving distance)
- $b(\chi)$ = linear bias
- $n(\chi)$ = number density of galaxies
- $D(\chi)$ = growth function

```julia
"""
    density_weight_function(chi::Float64, background::Background, 
                           survey_params::SurveyParams)::Float64

Weight function for density probe: f^den(χ) = H(χ)/c · b(χ) · n(χ) · D(χ)
"""
function density_weight_function(chi::Float64, background::Background, 
                                survey_params::SurveyParams)::Float64
    z = background.interp_chi_inv(chi)
    
    H_z = background.interp_H(z)  # 1/Mpc
    c_over_Mpc = 1.0 / 299792.458  # c in Mpc/Mpc (speed of light)
    
    bias_z = survey_params.get_bias(z)
    n_z = survey_params.get_number_density(z)
    D_z = background.interp_D(z)
    
    return (H_z * c_over_Mpc) * bias_z * n_z * D_z
end

struct SurveyParams
    # Number density: dN/dz
    n_z::Function
    
    # Bias function: b(z)
    b_z::Function
    
    # Could add: magnification bias, evolution bias, etc.
end

function SurveyParams(n_z_file::String, b_z_file::String)
    # Load from precomputed files
    n_z = load_interp(n_z_file)
    b_z = load_interp(b_z_file)
    return SurveyParams(n_z, b_z)
end
```

For **lensing** probe (cosmic shear):
$$f^{\text{lens}}(\chi) = \frac{3}{2} \Omega_m \left(\frac{H_0}{c}\right)^2 \int_\chi^\infty d\chi' \frac{a(\chi')}{E(z(\chi'))} \, n(\chi') \, \frac{\chi - \chi'}{\chi}$$

(More complex; lensing kernel involves integral over source distribution)

#### 3-ii. **Angular Window Function $W_\ell^A(k, \chi)$**

$$W_\ell^A(k, \chi) = f^A(\chi) \, j_\ell(k \chi)$$

```julia
"""
    window_function(ell::Int, k::Float64, chi::Float64,
                   f_weight::Function)::Float64

Angular window function: W_ℓ^A(k,χ) = f^A(χ) · j_ℓ(k·χ)
"""
function window_function(ell::Int, k::Float64, chi::Float64,
                        f_weight::Function)::Float64
    j_ell = spherical_bessel(ell, k * chi)
    return f_weight(chi) * j_ell
end
```

#### 3-iii. **Hankel Transform $\widetilde{W}_\ell^A(k, k_1)$**

The key computational step:

$$\widetilde{W}_\ell^A(k, k_1) = \int_0^\infty d\chi \, \chi^2 \, W_\ell^A(k, \chi) \, j_\ell(k_1 \chi)$$

$$= \int_0^\infty d\chi \, \chi^2 \, f^A(\chi) \, j_\ell(k\chi) \, j_\ell(k_1 \chi)$$

```julia
"""
    hankel_transform(ell::Int, k::Float64, k1::Float64,
                    z_arr::Vector{Float64}, f_weight::Function,
                    background::Background,
                    chi_interp::Interpolation)::Float64

Compute Hankel transform: W̃_ℓ(k, k₁)
Integrates over comoving distance χ ∈ [χ_min, χ_max]
Uses adaptive quadrature
"""
function hankel_transform(ell::Int, k::Float64, k1::Float64,
                         z_arr::Vector{Float64}, f_weight::Function,
                         background::Background,
                         chi_interp::Interpolation)::Float64
    
    # Get comoving distance range from redshift array
    chi_min = chi_interp(z_arr[1])
    chi_max = chi_interp(z_arr[end])
    
    # Integrand: χ² · f(χ) · j_ℓ(kχ) · j_ℓ(k₁χ)
    function integrand(chi::Float64)::Float64
        chi_sq = chi^2
        j_ell_k = spherical_bessel(ell, k * chi)
        j_ell_k1 = spherical_bessel(ell, k1 * chi)
        f_chi = f_weight(chi)
        
        return chi_sq * f_chi * j_ell_k * j_ell_k1
    end
    
    # Adaptive quadrature
    result, _ = quadgk(integrand, chi_min, chi_max, rtol=1e-5)
    
    return result
end

# Optimized version: vectorize over k₁ values
"""
    hankel_transform_grid(ell::Int, k::Float64, k1_arr::Vector{Float64},
                         z_arr::Vector{Float64}, f_weight::Function,
                         background::Background,
                         chi_interp::Interpolation)::Vector{Float64}

Compute W̃_ℓ(k, k₁) for array of k₁ values (vectorized)
"""
function hankel_transform_grid(ell::Int, k::Float64, k1_arr::Vector{Float64},
                              z_arr::Vector{Float64}, f_weight::Function,
                              background::Background,
                              chi_interp::Interpolation)::Vector{Float64}
    
    return [hankel_transform(ell, k, k1, z_arr, f_weight, background, chi_interp)
            for k1 in k1_arr]
end
```

#### 3-iv. **Cross-Hankel Transforms for Cross-Correlations**

For cross-correlations between two different probes (A and B):

$$\widetilde{W}_\ell^A(k, k_1) \times \widetilde{W}_\ell^B(k, k_2) = \int d\chi \chi^2 f^A(\chi) j_\ell(k\chi) j_\ell(k_1\chi) \times \int d\chi' (\chi')^2 f^B(\chi') j_\ell(k\chi') j_\ell(k_2\chi')$$

These are **separable** and can be precomputed independently.

**File location**: `src/sfb_kernels.jl` (NEW FILE)

---

### LAYER 4: SFB Coefficients (`sfb_decomposition.jl` - NEW FILE)

**Purpose**: Assemble the Hankel transforms into the final correlation function.

#### 4-i. **3D Correlation Function $S_\ell^{AB}(k_1, k_2)$**

From `equations.tex`, Eq.(ssfb_new):

$$S_\ell^{AB}(k_1, k_2) = N_\ell^{AB} \int dk \, k^2 \, P_{\text{lin}}^{AB}(k) \, \widetilde{W}_\ell^A(k, k_1) \, \widetilde{W}_\ell^B(k, k_2)$$

**Conceptual structure**:
1. Define grid in $(k, k_1, k_2)$ space
2. Precompute $\widetilde{W}_\ell^A$ and $\widetilde{W}_\ell^B$ on their respective $k$ grids
3. For each $(k_1, k_2)$ pair, integrate over $k$

```julia
"""
    struct SFBCorrelation
        ell::Int
        k_vals::Vector{Float64}              # Integration grid
        k1_vals::Vector{Float64}             # Output grid for k₁
        k2_vals::Vector{Float64}             # Output grid for k₂
        
        W_tilde_A::Matrix{Float64}           # W̃^A(k, k₁) [n_k × n_k1]
        W_tilde_B::Matrix{Float64}           # W̃^B(k, k₂) [n_k × n_k2]
        P_lin::Vector{Float64}               # P_lin(k)
    end
"""

"""
    compute_correlation_function(ell::Int,
                                k1::Float64, k2::Float64,
                                sfb_corr::SFBCorrelation,
                                cosmo::CosmoParams)::Float64

Compute S_ℓ^{AB}(k₁, k₂) via 1D quadrature over k
"""
function compute_correlation_function(ell::Int,
                                     k1::Float64, k2::Float64,
                                     sfb_corr::SFBCorrelation,
                                     cosmo::CosmoParams)::Float64
    
    # Interpolate W̃ at k₁, k₂
    W_tilde_A_at_k1 = interpolate(sfb_corr.k_vals, sfb_corr.W_tilde_A[:, :], k1)
    W_tilde_B_at_k2 = interpolate(sfb_corr.k_vals, sfb_corr.W_tilde_B[:, :], k2)
    
    # Integrand: k² · P_lin(k) · W̃^A(k,k₁) · W̃^B(k,k₂)
    function integrand(k::Float64)::Float64
        k_sq = k^2
        P_k = sfb_corr.P_lin(k)
        W_A = W_tilde_A_at_k1(k)
        W_B = W_tilde_B_at_k2(k)
        
        return k_sq * P_k * W_A * W_B
    end
    
    # Adaptive quadrature
    result, _ = quadgk(integrand, sfb_corr.k_vals[1], sfb_corr.k_vals[end], rtol=1e-5)
    
    # Normalization (from survey)
    N_AB = 1.0  # Depends on survey specifications
    
    return N_AB * result
end

"""
    compute_correlation_grid(ell::Int,
                            k1_arr::Vector{Float64}, k2_arr::Vector{Float64},
                            sfb_corr::SFBCorrelation,
                            cosmo::CosmoParams)::Matrix{Float64}

Compute S_ℓ^{AB} on full (k₁, k₂) grid (vectorized)
"""
function compute_correlation_grid(ell::Int,
                                 k1_arr::Vector{Float64}, k2_arr::Vector{Float64},
                                 sfb_corr::SFBCorrelation,
                                 cosmo::CosmoParams)::Matrix{Float64}
    
    n_k1 = length(k1_arr)
    n_k2 = length(k2_arr)
    S_ell = zeros(n_k1, n_k2)
    
    for (i, k1) in enumerate(k1_arr)
        for (j, k2) in enumerate(k2_arr)
            S_ell[i, j] = compute_correlation_function(ell, k1, k2, sfb_corr, cosmo)
        end
    end
    
    return S_ell
end
```

#### 4-ii. **Precomputation & Caching**

```julia
"""
    precompute_sfb_correlations(ell_arr::Vector{Int},
                               survey::SurveyParams,
                               background::Background,
                               P_lin::Function,
                               cosmo::CosmoParams;
                               n_k_integration::Int=100,
                               k_min::Float64=1e-3,
                               k_max::Float64=10.0)::Dict

Precompute all Hankel transforms and correlation functions for given ℓ values
"""
function precompute_sfb_correlations(ell_arr::Vector{Int},
                                    survey::SurveyParams,
                                    background::Background,
                                    P_lin::Function,
                                    cosmo::CosmoParams;
                                    n_k_integration::Int=100,
                                    k_min::Float64=1e-3,
                                    k_max::Float64=10.0)::Dict
    
    cache = Dict()
    
    for ell in ell_arr
        # Setup k grid
        k_vals = exp.(range(log(k_min), log(k_max), length=n_k_integration))
        
        # Compute Hankel transforms
        W_tilde_A = hankel_transform_grid(ell, k_vals, k_vals, ...)
        W_tilde_B = hankel_transform_grid(ell, k_vals, k_vals, ...)
        
        cache[ell] = SFBCorrelation(ell, k_vals, k_vals, k_vals,
                                    W_tilde_A, W_tilde_B, P_lin.(k_vals))
    end
    
    return cache
end
```

**File location**: `src/sfb_decomposition.jl` (NEW FILE)

---

### LAYER 5: Final Observables (`integrals.jl`, potentially expand)

**Purpose**: Integrate the 3D correlation function over redshift bins to get the final angular power spectrum $C_\ell^{AB}$.

#### 5-i. **From $S_\ell^{AB}(k_1, k_2)$ to $C_\ell^{AB}$**

The angular power spectrum is obtained by integrating $S_\ell$ over the three wavenumbers that appear (corresponding to three different redshift bins in 3x2pt):

$$C_\ell^{AB}(z_i, z_j) = \int \frac{dk}{2\pi^2} \int \frac{dk_1}{2\pi^2} \int \frac{dk_2}{2\pi^2} \, S_\ell^{AB}(k, k_1, k_2)$$

For practical computation, this becomes a **multidimensional integral** that requires careful handling.

#### 5-ii. **Angular Power Spectrum Function**

```julia
"""
    angular_power_spectrum(ell::Int,
                          z_i::Float64, z_j::Float64,
                          sfb_cache::Dict,
                          cosmo::CosmoParams)::Float64

Compute C_ℓ^{AB}(z_i, z_j) by integrating S_ℓ over wavenumbers
"""
function angular_power_spectrum(ell::Int,
                               z_i::Float64, z_j::Float64,
                               sfb_cache::Dict,
                               cosmo::CosmoParams)::Float64
    
    if !haskey(sfb_cache, ell)
        error("ℓ = $ell not in precomputed cache")
    end
    
    sfb_corr = sfb_cache[ell]
    
    # Triple integral over k, k₁, k₂
    function integrand(k_arr)
        k, k1, k2 = k_arr
        S_ell = compute_correlation_function(ell, k1, k2, sfb_corr, cosmo)
        # Factor from k³ (spherical volume element)
        return S_ell / (2π)^3
    end
    
    # Use high-dimensional quadrature (e.g., Cubature.jl)
    result, _ = hcubature(integrand, [k_min, k_min, k_min],
                                     [k_max, k_max, k_max],
                                     rtol=1e-4)
    
    return result
end
```

#### 5-iii. **Validation Against Limber Reference**

**Critical step**: Compare against precomputed Limber results to validate correctness.

```julia
"""
    validate_against_limber(ell_arr::Vector{Int},
                           z_bins::Vector{Tuple{Float64, Float64}},
                           C_ell_sfb::Matrix{Float64},
                           C_ell_limber::Matrix{Float64})::Dict

Compute residuals and convergence metrics
"""
function validate_against_limber(ell_arr::Vector{Int},
                                z_bins::Vector{Tuple{Float64, Float64}},
                                C_ell_sfb::Matrix{Float64},
                                C_ell_limber::Matrix{Float64})::Dict
    
    n_ell = length(ell_arr)
    n_bins = length(z_bins)
    
    relative_error = similar(C_ell_sfb)
    for i in 1:n_ell
        for j in 1:n_bins
            relative_error[i, j] = abs(C_ell_sfb[i, j] - C_ell_limber[i, j]) / 
                                    abs(C_ell_limber[i, j])
        end
    end
    
    return Dict(
        "max_rel_error" => maximum(relative_error),
        "mean_rel_error" => mean(relative_error),
        "relative_error_grid" => relative_error
    )
end
```

**File location**: `src/integrals.jl` (extend existing file)

---

## Mathematical Foundations

### Key Equations From `equations.tex`

#### Eq. 1: Generalized BLAST Formula
$$C_\ell^{AB} = N_\ell^{AB} \int d\chi \int d\chi' \int dk \, k^2 \, P^{AB}(k,\chi,\chi') \, j_\ell(k\chi) \, j_\ell(k\chi')$$

**Physical meaning**: 
- $C_\ell$ is the angular power spectrum (observable)
- $k$ is the Fourier wavenumber
- $\chi$ is the comoving distance (related to redshift)
- $j_\ell$ are spherical Bessel functions encoding the angular information

#### Eq. 2: Hankel Transform Strategy
$$\widetilde{W}_\ell^A(k, k_1) = \int_0^\infty d\chi \, \chi^2 \, W_\ell^A(k, \chi) \, j_\ell(k_1 \chi)$$

**Innovation**: This separates the integral, allowing efficient computation.

#### Eq. 3: 3D Power Spectrum Expansion
$$S_\ell^{AB}(k_1, k_2) = N_\ell^{AB} \int dk \, k^2 \, P_{\text{lin}}^{AB}(k) \, \widetilde{W}_\ell^A(k, k_1) \, \widetilde{W}_\ell^B(k, k_2)$$

**Advantage**: Now only 1D integration over $k$ (instead of 3D).

#### Eq. 4: Density Window Function
$$f^{\text{den}}(\chi) = \frac{H(\chi)}{c} \cdot b(\chi) \cdot n(\chi) \cdot D(\chi)$$

**Physical components**:
- $H(\chi)/c$: converts comoving distance to dynamical time
- $b(\chi)$: linear bias (how galaxies trace matter)
- $n(\chi)$: galaxy number density
- $D(\chi)$: growth factor (how perturbations grow)

---

## Code Structure & File Mapping

### Complete Module Organization

```
src/Blast.jl
├── cosmo.jl                    # [LAYER 0] Cosmological parameters
├── background.jl               # [LAYER 1] Background evolution: E(z), χ(z), D(z)
├── projected_matter.jl         # [LAYER 2A] Matter power spectra
├── spherical_bessel.jl         # [LAYER 2B] j_ℓ(x), zeros [NEW]
├── sfb_kernels.jl              # [LAYER 3] Window functions, Hankel transforms [NEW]
├── sfb_decomposition.jl        # [LAYER 4] S_ℓ(k₁,k₂) and SFB coefficients [NEW]
└── integrals.jl                # [LAYER 5] Final C_ℓ computation
```

### Import Hierarchy

```julia
# At module level (Blast.jl):
include("cosmo.jl")
include("background.jl")
include("projected_matter.jl")
include("spherical_bessel.jl")
include("sfb_kernels.jl")
include("sfb_decomposition.jl")
include("integrals.jl")
```

---

## Detailed Implementation Steps

### PHASE 1: Foundation (Weeks 1-2)

**Goal**: Set up core infrastructure that everything else depends on.

1. **Expand `cosmo.jl`**
   - [ ] Define `CosmoParams` struct with validation
   - [ ] Implement derived quantities: $H_0$, $\Omega_k$
   - [ ] Add docstrings and tests

2. **Expand `background.jl`**
   - [ ] Implement $E(z)$ function
   - [ ] Implement $\chi(z)$ with precomputation
   - [ ] Implement $D(z)$ from ODE solver
   - [ ] Create `Background` struct for caching
   - [ ] Validate against reference values

3. **Create test notebook**: `sfb_test_background.ipynb`
   - [ ] Plot $E(z)$, $\chi(z)$, $D(z)$ for standard cosmology
   - [ ] Compare $\chi(z)$ with precomputed data in `data/background/`

---

### PHASE 2: Input Data (Weeks 2-3)

**Goal**: Load and validate power spectrum and precomputed data.

1. **Expand `projected_matter.jl`**
   - [ ] Load $P_{\text{lin}}(k)$ from CLASS/CAMB or file
   - [ ] Implement Halofit or other nonlinear prescription
   - [ ] Implement $P^{3D}(k, \chi, \chi')$ with growth weighting
   - [ ] Vectorize for efficiency

2. **Create `SurveyParams` struct**
   - [ ] Load galaxy number density $n(z)$
   - [ ] Load bias $b(z)$
   - [ ] Add validation

3. **Test notebook**: `sfb_test_power_spectrum.ipynb`
   - [ ] Plot $P_{\text{lin}}(k)$ and $P_{\text{nl}}(k)$
   - [ ] Compare growth-weighted $P^{3D}$ at different $(z, z')$ pairs

---

### PHASE 3: Spherical Bessel Functions (Weeks 3-4)

**Goal**: Implement efficient spherical Bessel function computation.

1. **Create `spherical_bessel.jl`**
   - [ ] Implement $j_\ell(x)$ using SpecialFunctions.jl
   - [ ] Implement derivatives $j'_\ell(x)$
   - [ ] Implement zero-finding for first $n$ zeros
   - [ ] Add caching for frequently used values
   - [ ] Vectorize all functions

2. **Optimize and benchmark**
   - [ ] Profile compute time for typical usage patterns
   - [ ] Compare against reference implementations (GSL, scipy)

3. **Test notebook**: `sfb_test_bessel.ipynb`
   - [ ] Plot $j_\ell(x)$ for several $\ell$ values
   - [ ] Verify zeros against reference values
   - [ ] Benchmark and compare with/without caching

---

### PHASE 4: Window Functions & Hankel Transforms (Weeks 4-6)

**Goal**: Implement the core Hankel transform machinery.

1. **Create `sfb_kernels.jl`**
   - [ ] Implement density weight function $f^{\text{den}}(\chi)$
   - [ ] Implement lensing weight function (if included)
   - [ ] Implement angular window function $W_\ell^A(k, \chi)$
   - [ ] Implement Hankel transform $\widetilde{W}_\ell^A(k, k_1)$
   - [ ] Vectorize Hankel transform over $k_1$ grid

2. **Validation**
   - [ ] Check Hankel transform against analytical limits (if available)
   - [ ] Verify integrals converge properly
   - [ ] Test with simple weight functions (e.g., constant)

3. **Optimization**
   - [ ] Profile quadrature time
   - [ ] Tune tolerances and grid spacing
   - [ ] Consider caching strategies

4. **Test notebook**: `sfb_test_hankel.ipynb`
   - [ ] Plot $\widetilde{W}_\ell^A(k, k_1)$ for several multipoles
   - [ ] Verify convergence of integral
   - [ ] Compare Hankel transforms for different window functions

---

### PHASE 5: 3D Correlation Function (Weeks 6-7)

**Goal**: Assemble correlation functions from precomputed Hankel transforms.

1. **Create `sfb_decomposition.jl`**
   - [ ] Define `SFBCorrelation` struct
   - [ ] Implement $S_\ell^{AB}(k_1, k_2)$ computation
   - [ ] Implement grid computation (vectorized over $(k_1, k_2)$ pairs)
   - [ ] Implement precomputation and caching

2. **Validation Against Limber**
   - [ ] Compute $S_\ell$ for test cases
   - [ ] Load precomputed Limber reference from `data/Limber/`
   - [ ] Compare deviations (should be small for small scales)
   - [ ] Document convergence vs number of SFB modes

3. **Performance analysis**
   - [ ] Measure wall time and memory usage
   - [ ] Compare scaling with $\ell$, $k$-grid size
   - [ ] Identify bottlenecks

4. **Test notebook**: `sfb_test_correlation.ipynb`
   - [ ] Plot $S_\ell^{AB}(k_1, k_2)$ as 2D heatmaps
   - [ ] Compare with Limber reference
   - [ ] Show convergence plots

---

### PHASE 6: Final Angular Power Spectrum (Weeks 7-8)

**Goal**: Complete the pipeline by computing final observables.

1. **Extend `integrals.jl`**
   - [ ] Implement $C_\ell^{AB}$ computation from $S_\ell$
   - [ ] Handle multidimensional integration (high-dimensional quadrature)
   - [ ] Implement validation against Limber reference

2. **Validation**
   - [ ] Compute $C_\ell$ for all redshift bin pairs
   - [ ] Compare with `data/Limber/Cl_*.npy` precomputed data
   - [ ] Measure residuals and convergence

3. **Documentation**
   - [ ] Write README for new modules
   - [ ] Document equations-to-code mapping
   - [ ] Provide usage examples

4. **Test notebook**: `sfb_development.ipynb` (final)
   - [ ] End-to-end pipeline test
   - [ ] Generate comparison plots with Limber
   - [ ] Demonstrate 3x2pt analysis results

---

### PHASE 7: Optimization & Scaling (Weeks 8+)

**Goal**: Make code production-ready with high performance.

1. **Optimize critical paths**
   - [ ] Use `@turbo` for inner loops (LoopVectorization.jl)
   - [ ] Use `@tullio` for tensor contractions (Tullio.jl)
   - [ ] Profile with BenchmarkTools
   - [ ] Compare scaling: current vs BLAST vs other codes

2. **Add advanced features**
   - [ ] Support for different nonlinear spectra (HMCode, EuclidEmulator, etc.)
   - [ ] Support for magnification bias and other effects
   - [ ] Batch computation for multiple cosmologies
   - [ ] Error propagation framework

3. **Production validation**
   - [ ] Test suite for all modules
   - [ ] Regression tests against Limber
   - [ ] Documentation of algorithms

---

## Summary: Execution Roadmap

| Phase | Duration | Key Files | Milestones |
|-------|----------|-----------|-----------|
| 1 | Weeks 1-2 | cosmo.jl, background.jl | Background evolution precomputed & validated |
| 2 | Weeks 2-3 | projected_matter.jl | Power spectrum loaded & cached |
| 3 | Weeks 3-4 | spherical_bessel.jl | Bessel functions optimized & tested |
| 4 | Weeks 4-6 | sfb_kernels.jl | Hankel transforms working & benchmarked |
| 5 | Weeks 6-7 | sfb_decomposition.jl | S_ℓ(k₁,k₂) validated vs Limber |
| 6 | Weeks 7-8 | integrals.jl | C_ℓ^AB computed end-to-end |
| 7 | Weeks 8+ | all files | Full optimization & production deployment |

---

## Next Steps

1. **Start with Phase 1**: Expand `background.jl` to compute and cache $\chi(z)$ and $D(z)$
2. **Create validation notebook** to confirm outputs match precomputed data
3. **Incrementally add layers**, testing each against references
4. **Document equations-to-code** mapping at each step (included above)

This guide provides a complete blueprint for translating the math in `equations.tex` into efficient, production-ready Julia code following the BLAST framework.

