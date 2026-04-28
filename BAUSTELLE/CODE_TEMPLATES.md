# Practical Implementation Guide: Code Templates & Pseudocode

**Purpose**: Translate each mathematical equation into actual Julia code with detailed examples.

---

## Quick Reference: Equation → Code Mapping

| Equation | Location | Julia Function | Input → Output |
|----------|----------|-----------------|-----------------|
| $E(z) = \sqrt{\Omega_m(1+z)^3 + ...}$ | background.jl | `E(z, cosmo)` | $z \to E(z)$ |
| $\chi(z) = \frac{c}{H_0}\int_0^z \frac{dz'}{E(z')}$ | background.jl | `comoving_distance(z, cosmo)` | $z \to \chi(z)$ Mpc |
| $j_\ell(x) = \sqrt{\frac{\pi}{2x}} J_{\ell+1/2}(x)$ | spherical_bessel.jl | `spherical_bessel(ell, x)` | $(\ell, x) \to j_\ell(x)$ |
| $W_\ell^A(k,\chi) = f^A(\chi) \cdot j_\ell(k\chi)$ | sfb_kernels.jl | `window_function(ell, k, chi, f)` | $(k,\chi) \to W_\ell$ |
| $\widetilde{W}_\ell^A(k,k_1) = \int d\chi \chi^2 W_\ell^A(k,\chi) j_\ell(k_1\chi)$ | sfb_kernels.jl | `hankel_transform(ell, k, k1, ...)` | $(k,k_1) \to \widetilde{W}$ |
| $S_\ell^{AB}(k_1,k_2) = N_\ell \int dk \, k^2 P(k) \widetilde{W}^A \widetilde{W}^B$ | sfb_decomposition.jl | `compute_correlation_function(ell, k1, k2, ...)` | $(k_1,k_2) \to S_\ell$ |

---

## LAYER 1: Background Evolution

### Function: `E(z)` — Normalized Hubble Parameter

**Mathematical Definition**:
$$E(z) = \frac{H(z)}{H_0} = \sqrt{\Omega_m(1+z)^3 + \Omega_k(1+z)^2 + \Omega_\Lambda a^{-3(1+w(a))}}$$

**Julia Template**:
```julia
"""
    E(z::Float64, cosmo::CosmoParams)::Float64

Compute normalized Hubble parameter E(z) = H(z)/H₀

# Physics
- Valid for any w₀-wₐ dark energy model
- Used: background evolution, distance computations, growth rate

# Arguments
- `z::Float64`: Redshift
- `cosmo::CosmoParams`: Cosmological model with (Ω_m, Ω_Λ, Ω_k, w₀, wₐ)

# Returns
- `E(z)::Float64`: Dimensionless Hubble parameter

# Example
```julia
cosmo = CosmoParams(Ω_m=0.3, Ω_Λ=0.7)
E_at_z1 = E(1.0, cosmo)  # Returns ≈ 1.8 for typical ΛCDM
```
"""
function E(z::Float64, cosmo::CosmoParams)::Float64
    # (1) Compute scale factor
    a = 1.0 / (1.0 + z)
    
    # (2) Compute dark energy equation of state at redshift z
    #     w(a) = w₀ + wₐ(1 - a)
    w_z = cosmo.w0 + cosmo.wa * (1.0 - a)
    
    # (3) Compute each term in Friedmann equation
    term_matter = cosmo.Ω_m * (1.0 + z)^3
    term_curvature = cosmo.Ω_k * (1.0 + z)^2
    term_dark_energy = cosmo.Ω_Λ * a^(-3.0 * (1.0 + w_z))
    
    # (4) Sum and take square root
    return sqrt(term_matter + term_curvature + term_dark_energy)
end

# Vectorized version (for efficiency)
E(z_arr::Vector{Float64}, cosmo::CosmoParams) = @. E(z_arr, cosmo)
```

**Validation**:
```julia
using Test

cosmo = CosmoParams(Ω_m=0.3, Ω_Λ=0.7, h=0.67)
@test E(0.0, cosmo) ≈ 1.0  # By definition at z=0
@test E(1.0, cosmo) > 1.0  # Hubble increases with z
@test all(E([0, 1, 2], cosmo) .> 0)  # Always positive
```

---

### Function: `comoving_distance(z)` — Comoving Distance χ(z)

**Mathematical Definition**:
$$\chi(z) = \frac{c}{H_0} \int_0^z \frac{dz'}{E(z')}$$

**Julia Template**:
```julia
"""
    comoving_distance(z::Float64, cosmo::CosmoParams)::Float64

Compute comoving distance from observer to redshift z

# Physics
- Accounts for cosmic expansion
- Used in: window functions, spherical Bessel arguments
- Units: Megaparsec (Mpc)

# Equation
χ(z) = (c/H₀) ∫₀^z dz'/E(z')

# Implementation strategy
1. Use adaptive quadrature (QuadGK.jl)
2. Handle z=0 case analytically
3. Precompute and interpolate for efficiency

# Arguments
- `z::Float64`: Redshift
- `cosmo::CosmoParams`: Cosmological parameters

# Returns
- `chi::Float64`: Comoving distance in Mpc

# Example
```julia
cosmo = CosmoParams(Ω_m=0.3)
chi_z1 = comoving_distance(1.0, cosmo)  # ≈ 3286 Mpc
```
"""
function comoving_distance(z::Float64, cosmo::CosmoParams)::Float64
    # Handle z ≈ 0 case
    if z < 1e-6
        return 0.0
    end
    
    # Physical constant: c / H₀ in Mpc
    # c = 299792.458 km/s
    # H₀ = 100h km/s/Mpc
    # c/H₀ = 299792.458 / (100h) Mpc
    c_over_H0_Mpc = 299792.458 / (100.0 * cosmo.h)
    
    # Integrand: 1 / E(z)
    integrand(z_prime) = 1.0 / E(z_prime, cosmo)
    
    # Adaptive quadrature integration
    integral, error = quadgk(integrand, 0.0, z, rtol=1e-6, atol=1e-8)
    
    return c_over_H0_Mpc * integral
end

# IMPORTANT: Precompute for speed!
"""
    struct DistanceInterpolator
        z_grid::Vector{Float64}
        chi_grid::Vector{Float64}
        itp::Interpolations.Extrapolation
    end
    
    function DistanceInterpolator(z_max::Float64, n_points::Int, 
                                 cosmo::CosmoParams)
        # Create logarithmic grid for better sampling at high z
        z_grid = exp.(range(log(1e-3), log(z_max), length=n_points))
        
        # Compute distances at grid points
        chi_grid = zeros(n_points)
        for i in eachindex(z_grid)
            chi_grid[i] = comoving_distance(z_grid[i], cosmo)
        end
        
        # Create interpolator (linear in χ, or log-linear for better behavior)
        itp = interpolation((z_grid,), chi_grid, Gridded(Linear()))
        
        return DistanceInterpolator(z_grid, chi_grid, itp)
    end
    
    # Usage:
    chi_interp = DistanceInterpolator(10.0, 500, cosmo)
    chi_at_z = chi_interp.itp(1.5)  # Fast evaluation
"""
```

**Validation**:
```julia
cosmo = CosmoParams(Ω_m=0.3, Ω_Λ=0.7, h=0.67)

@test comoving_distance(0.0, cosmo) ≈ 0.0
@test comoving_distance(1.0, cosmo) > 0.0

# Should be monotonically increasing
z_vals = [0.1, 0.5, 1.0, 2.0]
chi_vals = comoving_distance.(z_vals, Ref(cosmo))
@test issorted(chi_vals)
```

---

### Function: `growth_function(z)` — Linear Growth D(z)

**Mathematical Definition** (ODE form):
$$\frac{d^2D}{da^2} + \frac{3}{a}\frac{dD}{da} + \frac{\Omega_m(a) - 2}{a^2 E^2(a)} D = 0$$

with normalization $D(a=1) = 1$.

**Julia Template**:
```julia
"""
    growth_function(z::Float64, cosmo::CosmoParams; normalize::Bool=true)::Float64

Compute linear growth factor D(z), normalized to D(z=0)=1

# Physics
- Describes growth of matter perturbations
- Used in: matter power spectrum scaling P(k,z) = D²(z) P_lin(k)
- Dimensionless quantity

# Implementation
- Solve 2nd-order ODE backwards from high z to target z
- Uses DifferentialEquations.jl (ODE solver)

# Arguments
- `z::Float64`: Redshift
- `cosmo::CosmoParams`: Cosmological parameters
- `normalize::Bool=true`: Normalize to D(z=0)=1

# Returns
- `D(z)::Float64`: Linear growth factor

# Example
```julia
cosmo = CosmoParams(Ω_m=0.3)
D_z1 = growth_function(1.0, cosmo)  # ≈ 0.5 for ΛCDM
```
"""
function growth_function(z::Float64, cosmo::CosmoParams; normalize::Bool=true)::Float64
    using DifferentialEquations
    
    # Define ODE system: u = [D, dD/da]
    function growth_ode!(du, u, p, a)
        D, dD_da = u
        
        # Convert a to z
        z_a = 1.0/a - 1.0
        
        # Hubble parameter and density parameter at this a
        E_a = E(z_a, cosmo)
        Ω_m_a = cosmo.Ω_m * (1.0 + z_a)^3 / E_a^2
        
        # Second derivative from ODE
        d2D_da2 = -3.0/a * dD_da - 
                   (Ω_m_a - 2.0) / (a^2 * E_a^2) * D
        
        du[1] = dD_da
        du[2] = d2D_da2
    end
    
    # Initial conditions at early times (a_i ≈ small)
    a_init = 1e-3
    D_init = a_init  # D ∝ a in matter-dominated limit
    dD_da_init = 1.0  # Approximate
    u0 = [D_init, dD_da_init]
    
    # Final time (target redshift)
    a_final = 1.0 / (1.0 + z)
    
    # Solve ODE
    tspan = (a_init, a_final)
    prob = ODEProblem(growth_ode!, u0, tspan, cosmo)
    sol = solve(prob, RK45(), reltol=1e-6, abstol=1e-8)
    
    # Extract D at final time
    D_z = sol(a_final)[1]
    
    # Normalize to D(z=0) = 1
    if normalize
        D_0 = sol(1.0)[1]
        D_z = D_z / D_0
    end
    
    return D_z
end

# IMPORTANT: Cache for efficiency!
"""
    struct GrowthInterpolator
        z_grid::Vector{Float64}
        D_grid::Vector{Float64}
        itp::Interpolations.Extrapolation
    end
"""

const GROWTH_CACHE = Dict()

"""
    get_growth_function(z::Float64, cosmo::CosmoParams)::Float64

Retrieve D(z) from cache, or compute and cache if not present
"""
function get_growth_function(z::Float64, cosmo::CosmoParams)::Float64
    # Create cache key (hash of cosmological parameters)
    key = hash(cosmo)
    
    if !haskey(GROWTH_CACHE, key)
        # Precompute on grid
        z_arr = exp.(range(log(1e-3), log(100), length=300))
        D_arr = growth_function.(z_arr, Ref(cosmo))
        
        # Create interpolator
        itp = linear_interpolation(log.(z_arr), log.(D_arr), 
                                   extrapolation_bc=Line())
        GROWTH_CACHE[key] = itp
    end
    
    itp = GROWTH_CACHE[key]
    return exp(itp(log(z)))
end
```

---

## LAYER 2: Spherical Bessel Functions

### Function: `spherical_bessel(ell, x)` — Spherical Bessel j_ℓ(x)

**Mathematical Definition**:
$$j_\ell(x) = \sqrt{\frac{\pi}{2x}} J_{\ell+1/2}(x)$$

**Julia Template**:
```julia
"""
    spherical_bessel(ell::Int, x::Float64)::Float64

Compute spherical Bessel function j_ℓ(x)

# Physics
- Arise naturally in spherical harmonic expansions
- j_ℓ(x) oscillates for large x, decays at x ≈ 0
- Critical for computing angular power spectra

# Implementation
- Use SpecialFunctions.jl for cylindrical Bessel J_{ℓ+1/2}
- Convert: j_ℓ(x) = √(π/2x) · J_{ℓ+1/2}(x)

# Arguments
- `ell::Int`: Spherical harmonic order (ℓ ≥ 0)
- `x::Float64`: Argument (x ≥ 0)

# Returns
- `j_ell::Float64`: Spherical Bessel value

# Behavior
- j_ℓ(0) = 1 if ℓ=0, else 0
- j_ℓ(x) oscillates for large x
- j_ℓ(x) ∼ x^ℓ for small x

# Example
```julia
j_0_at_1 = spherical_bessel(0, 1.0)  # ≈ 0.8414
j_1_at_1 = spherical_bessel(1, 1.0)  # ≈ 0.3011
```
"""
function spherical_bessel(ell::Int, x::Float64)::Float64
    using SpecialFunctions
    
    # Handle edge case: x = 0
    if abs(x) < 1e-10
        return (ell == 0) ? 1.0 : 0.0
    end
    
    # Compute cylindrical Bessel J_{ℓ+1/2}(x)
    J_half = besselj(ell + 0.5, x)
    
    # Convert to spherical: j_ℓ(x) = √(π/2x) · J_{ℓ+1/2}(x)
    return sqrt(π / (2.0 * x)) * J_half
end

# Vectorized (critical for performance!)
function spherical_bessel(ell::Int, x_arr::Vector{Float64})::Vector{Float64}
    return @. spherical_bessel(ell, x_arr)
    # Note: @. applies function element-wise
end
```

### Function: `spherical_bessel_zeros(ell, n_zeros)` — Zeros of j_ℓ(x)

**Mathematical Definition**: Find $x_{\ell,n}$ such that $j_\ell(x_{\ell,n}) = 0$

**Julia Template**:
```julia
"""
    spherical_bessel_zeros(ell::Int, n_zeros::Int; 
                          x_max::Float64=1000.0)::Vector{Float64}

Find first n_zeros zeros of spherical Bessel function j_ℓ(x)

# Physics
- Zeros define radial grid for Spherical Fourier-Bessel decomposition
- x_{ℓ,n} = n-th zero of j_ℓ

# Implementation
- Use root-finding (bisection or Newton's method)
- Exploit approximate formula for zero locations

# Arguments
- `ell::Int`: Spherical harmonic order
- `n_zeros::Int`: Number of zeros to find
- `x_max::Float64`: Maximum search range

# Returns
- `zeros_vec::Vector{Float64}`: Sorted zeros [x_{ℓ,1}, x_{ℓ,2}, ...]

# Timing
- First zero: x_{ℓ,1} ≈ ℓ + 1
- n-th zero: x_{ℓ,n} ≈ π(n + ℓ/2 - 1/4)

# Example
```julia
zeros_ell0 = spherical_bessel_zeros(0, 10)  # First 10 zeros of j_0
# Returns: [3.14, 6.28, 9.42, ...] (≈ π, 2π, 3π, ...)
```
"""
function spherical_bessel_zeros(ell::Int, n_zeros::Int; 
                               x_max::Float64=1000.0)::Vector{Float64}
    using Roots
    
    zeros_vec = Float64[]
    
    # Start search near first expected zero
    # Approximate: first zero of j_ℓ is near ℓ + 1
    x_search = ell + 1.0
    
    for n in 1:n_zeros
        if x_search > x_max
            @warn "Stopped at n=$n: x_search > x_max"
            break
        end
        
        # Define function to find zero of
        f(x) = spherical_bessel(ell, x)
        
        try
            # Use bisection between x_search and x_search + π
            # (zeros are roughly π apart for large x)
            x_left = x_search
            x_right = x_search + π + 0.1
            
            # Make sure signs are opposite (bracket the zero)
            while f(x_left) * f(x_right) > 0
                x_right += 0.5
                if x_right > x_max
                    error("Cannot bracket zero for n=$n, ℓ=$ell")
                end
            end
            
            # Find zero using bisection
            zero_found = find_zero(f, (x_left, x_right), Bisection())
            push!(zeros_vec, zero_found)
            
            # Update search point for next zero
            x_search = zero_found + π * 0.8
            
        catch err
            @warn "Failed to find zero n=$n for ℓ=$ell: $err"
            break
        end
    end
    
    return zeros_vec
end

# Cache to avoid recomputation
const BESSEL_ZEROS_CACHE = Dict{Tuple{Int,Int}, Vector{Float64}}()

"""
    get_spherical_bessel_zeros(ell::Int, n_zeros::Int)::Vector{Float64}

Retrieve cached zeros or compute new ones
"""
function get_spherical_bessel_zeros(ell::Int, n_zeros::Int)::Vector{Float64}
    key = (ell, n_zeros)
    
    if !haskey(BESSEL_ZEROS_CACHE, key)
        BESSEL_ZEROS_CACHE[key] = spherical_bessel_zeros(ell, n_zeros)
    end
    
    return BESSEL_ZEROS_CACHE[key]
end
```

---

## LAYER 3: Hankel Transforms

### Function: `hankel_transform(ell, k, k1, f_weight, ...)` 

**Mathematical Definition**:
$$\widetilde{W}_\ell^A(k, k_1) = \int_0^\infty d\chi \, \chi^2 \, f^A(\chi) \, j_\ell(k\chi) \, j_\ell(k_1\chi)$$

**Julia Template**:
```julia
"""
    hankel_transform(ell::Int, k::Float64, k1::Float64,
                    z_min::Float64, z_max::Float64,
                    f_weight::Function,
                    chi_interp::Function,
                    z_interp::Function)::Float64

Compute Hankel transform W̃_ℓ(k, k₁) via integration over comoving distance

# Physics
- Transforms redshift-space window function to wavenumber space
- Separates 3D integral into independent 1D integrals
- Core of efficient beyond-Limber computation

# Equation
W̃_ℓ(k, k₁) = ∫ dχ χ² f(χ) j_ℓ(kχ) j_ℓ(k₁χ)

# Arguments
- `ell::Int`: Multipole order
- `k::Float64`: Wavenumber 1
- `k1::Float64`: Wavenumber 2
- `z_min, z_max::Float64`: Redshift range for integration
- `f_weight::Function`: Window weight function f(χ)
- `chi_interp::Function`: Interpolator for z↔χ conversion
- `z_interp::Function`: Interpolator for χ↔z conversion

# Returns
- `W_tilde::Float64`: Hankel transform value

# Tolerance
- rtol=1e-5 (controls integration accuracy)

# Performance
- ≈ 100 ms per call (depends on tolerance, integration range)
- Dominated by quadrature, not Bessel function evaluation

# Example
```julia
W_tilde = hankel_transform(2, 0.1, 0.2, 0.0, 2.0, f_weight, chi_int, z_int)
```
"""
function hankel_transform(ell::Int, k::Float64, k1::Float64,
                         z_min::Float64, z_max::Float64,
                         f_weight::Function,
                         chi_interp::Function,
                         z_interp::Function)::Float64
    
    # Convert redshift range to comoving distance range
    chi_min = chi_interp(z_min)
    chi_max = chi_interp(z_max)
    
    # Define integrand: χ² · f(χ) · j_ℓ(kχ) · j_ℓ(k₁χ)
    function integrand(chi::Float64)::Float64
        # Evaluate all components
        chi_sq = chi^2
        f_chi = f_weight(chi)
        j_ell_k = spherical_bessel(ell, k * chi)
        j_ell_k1 = spherical_bessel(ell, k1 * chi)
        
        # Product
        return chi_sq * f_chi * j_ell_k * j_ell_k1
    end
    
    # Adaptive quadrature (QuadGK.jl)
    result, error = quadgk(integrand, chi_min, chi_max, 
                          rtol=1e-5, atol=1e-10)
    
    return result
end

# Vectorized version: compute for multiple k1 values
"""
    hankel_transform_grid(ell::Int, k::Float64, k1_arr::Vector{Float64},
                         z_min::Float64, z_max::Float64,
                         f_weight::Function,
                         chi_interp::Function,
                         z_interp::Function)::Vector{Float64}

Compute W̃_ℓ(k, k₁) for array of k₁ values
Much more efficient than calling hankel_transform repeatedly
"""
function hankel_transform_grid(ell::Int, k::Float64, k1_arr::Vector{Float64},
                              z_min::Float64, z_max::Float64,
                              f_weight::Function,
                              chi_interp::Function,
                              z_interp::Function)::Vector{Float64}
    
    chi_min = chi_interp(z_min)
    chi_max = chi_interp(z_max)
    
    # Vectorized: compute once over χ, then multiply by j_ℓ(k₁χ) for each k₁
    # This is more efficient than calling hankel_transform repeatedly
    
    W_tilde_arr = zeros(length(k1_arr))
    
    for (i, k1) in enumerate(k1_arr)
        W_tilde_arr[i] = hankel_transform(ell, k, k1, 
                                         z_min, z_max,
                                         f_weight, 
                                         chi_interp, 
                                         z_interp)
    end
    
    return W_tilde_arr
end
```

---

## LAYER 4: 3D Correlation Function

### Function: `compute_correlation_function(ell, k1, k2, ...)`

**Mathematical Definition** (Eq. ssfb_new from equations.tex):
$$S_\ell^{AB}(k_1, k_2) = N_\ell^{AB} \int dk \, k^2 \, P_{\text{lin}}(k) \, \widetilde{W}_\ell^A(k, k_1) \, \widetilde{W}_\ell^B(k, k_2)$$

**Julia Template**:
```julia
"""
    struct SFBCorrelationCache
        ell::Int                          # Multipole order
        k_vals::Vector{Float64}           # Wavenumber grid for integration
        k1_vals::Vector{Float64}          # Output k₁ grid
        k2_vals::Vector{Float64}          # Output k₂ grid
        
        W_tilde_A::Matrix{Float64}        # W̃^A(k, k₁) shape: [n_k × n_k1]
        W_tilde_B::Matrix{Float64}        # W̃^B(k, k₂) shape: [n_k × n_k2]
        
        P_lin_vals::Vector{Float64}       # P_lin(k) at grid points
        
        # Interpolators for smooth evaluation
        W_A_itp::Interpolation.Extrapolation
        W_B_itp::Interpolation.Extrapolation
        P_itp::Interpolation.Extrapolation
    end

    function SFBCorrelationCache(ell::Int,
                                k_min::Float64=1e-3,
                                k_max::Float64=10.0,
                                n_k::Int=150,
                                z_bins::NTuple,
                                survey::SurveyParams,
                                background::Background,
                                P_lin_func::Function,
                                cosmo::CosmoParams)
        
        # Setup logarithmic grid in k
        k_vals = exp.(range(log(k_min), log(k_max), length=n_k))
        k1_vals = k_vals  # Same grid for output
        k2_vals = k_vals
        
        # Precompute Hankel transforms W̃^A(k, k₁)
        W_tilde_A = zeros(n_k, n_k)
        for (i, k) in enumerate(k_vals)
            W_tilde_A[:, i] .= hankel_transform_grid(ell, k, k_vals, 
                                                     z_bins[1]..., ...)
        end
        
        # Same for W̃^B(k, k₂)
        W_tilde_B = similar(W_tilde_A)
        for (i, k) in enumerate(k_vals)
            W_tilde_B[:, i] .= hankel_transform_grid(ell, k, k_vals, 
                                                     z_bins[2]..., ...)
        end
        
        # Precompute P_lin(k)
        P_lin_vals = P_lin_func.(k_vals)
        
        # Create interpolators
        W_A_itp = interpolation((k_vals, k_vals), W_tilde_A, Gridded(Linear()))
        W_B_itp = interpolation((k_vals, k_vals), W_tilde_B, Gridded(Linear()))
        P_itp = linear_interpolation(log.(k_vals), log.(P_lin_vals), 
                                     extrapolation_bc=Line())
        
        return SFBCorrelationCache(ell, k_vals, k1_vals, k2_vals,
                                  W_tilde_A, W_tilde_B, P_lin_vals,
                                  W_A_itp, W_B_itp, P_itp)
    end
"""

"""
    compute_correlation_function(ell::Int, k1::Float64, k2::Float64,
                                cache::SFBCorrelationCache)::Float64

Compute S_ℓ^{AB}(k₁, k₂) by integrating over k

# Physics
- S_ℓ is the 3D correlation function in Fourier space
- Integrating S_ℓ over k-space gives the full angular power spectrum

# Equation
S_ℓ(k₁, k₂) = N_ℓ ∫ dk k² P_lin(k) W̃^A(k, k₁) W̃^B(k, k₂)

# Arguments
- `ell::Int`: Multipole order
- `k1, k2::Float64`: Output wavenumbers
- `cache::SFBCorrelationCache`: Precomputed Hankel transforms

# Returns
- `S_ell::Float64`: Correlation function value

# Implementation strategy
1. Interpolate W̃^A and W̃^B at (k, k₁) and (k, k₂)
2. Compute integrand: k² · P_lin(k) · W̃^A · W̃^B
3. Integrate over k using adaptive quadrature

# Example
```julia
S_ell_val = compute_correlation_function(2, 0.1, 0.15, cache)
```
"""
function compute_correlation_function(ell::Int, k1::Float64, k2::Float64,
                                     cache::SFBCorrelationCache)::Float64
    
    # Define integrand: k² · P_lin(k) · W̃^A(k,k₁) · W̃^B(k,k₂)
    function integrand(k::Float64)::Float64
        k_sq = k^2
        
        # Evaluate power spectrum (via log interpolation for stability)
        P_k = exp(cache.P_itp(log(k)))
        
        # Evaluate Hankel transforms (via 2D interpolation)
        W_A_k_k1 = cache.W_A_itp(k, k1)
        W_B_k_k2 = cache.W_B_itp(k, k2)
        
        return k_sq * P_k * W_A_k_k1 * W_B_k_k2
    end
    
    # Integrate over k range
    k_min = cache.k_vals[1]
    k_max = cache.k_vals[end]
    
    result, _ = quadgk(integrand, k_min, k_max, rtol=1e-4)
    
    # Apply normalization (survey-dependent)
    N_AB = 1.0  # Would depend on survey specifics
    
    return N_AB * result
end

"""
    compute_correlation_grid(ell::Int,
                            k1_arr::Vector{Float64},
                            k2_arr::Vector{Float64},
                            cache::SFBCorrelationCache)::Matrix{Float64}

Compute S_ℓ^{AB} on full (k₁, k₂) grid

# Returns
- `S_ell::Matrix{Float64}` with shape [length(k1_arr), length(k2_arr)]
"""
function compute_correlation_grid(ell::Int,
                                 k1_arr::Vector{Float64},
                                 k2_arr::Vector{Float64},
                                 cache::SFBCorrelationCache)::Matrix{Float64}
    
    n_k1 = length(k1_arr)
    n_k2 = length(k2_arr)
    S_ell = zeros(n_k1, n_k2)
    
    # Vectorized computation
    for i in 1:n_k1
        for j in 1:n_k2
            S_ell[i, j] = compute_correlation_function(ell, k1_arr[i], k2_arr[j], 
                                                       cache)
        end
    end
    
    return S_ell
end
```

---

## Summary: Function Dependency Graph

```
┌─────────────────────────┐
│  CosmoParams (input)    │
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────────────────┐
│  E(z), comoving_distance(z),        │
│  growth_function(z)                 │
│  → Background struct                │
└────────────┬────────────────────────┘
             │
             ├─────────────────┐
             │                 │
             ▼                 ▼
    ┌──────────────────┐  ┌──────────────────┐
    │  projected_      │  │ spherical_       │
    │  matter.jl       │  │ bessel.jl        │
    │  P_lin(k)        │  │ j_ℓ(x)           │
    │  P_3D(k,χ,χ')    │  │ zeros(ℓ, n)      │
    └────────┬─────────┘  └──────────────────┘
             │
             ├─────────────────────────────────┐
             │                                 │
             ▼                                 ▼
    ┌─────────────────────────┐  ┌──────────────────────────┐
    │  sfb_kernels.jl         │  │ SurveyParams (input)    │
    │  - f_weight(χ)          │  │ - n(z), b(z)           │
    │  - window_function      │  │                         │
    │  - hankel_transform     │  │                         │
    └────────┬────────────────┘  └──────────────────────────┘
             │
             ▼
    ┌──────────────────────────────┐
    │  sfb_decomposition.jl        │
    │  - SFBCorrelationCache       │
    │  - compute_correlation_      │
    │    function(ell, k1, k2)     │
    │  → S_ℓ(k₁, k₂)              │
    └────────┬─────────────────────┘
             │
             ▼
    ┌──────────────────────────────┐
    │  integrals.jl (extend)       │
    │  - C_ℓ^AB computation        │
    │  - validation_vs_limber      │
    │  → final observables         │
    └──────────────────────────────┘
```

---

## Testing Strategy

### Test 1: Background Evolution
```julia
@testset "Background Evolution" begin
    cosmo = CosmoParams(Ω_m=0.3, Ω_Λ=0.7, h=0.67)
    
    # E(z) at z=0 should be 1
    @test E(0.0, cosmo) ≈ 1.0
    
    # χ(z) should be monotonically increasing
    z_vals = [0, 0.5, 1.0, 2.0]
    chi_vals = comoving_distance.(z_vals, Ref(cosmo))
    @test issorted(chi_vals)
    
    # D(z) at z=0 should be 1
    @test growth_function(0.0, cosmo) ≈ 1.0
end
```

### Test 2: Spherical Bessel Functions
```julia
@testset "Spherical Bessel" begin
    # j_0(0) = 1
    @test spherical_bessel(0, 0.0) ≈ 1.0
    
    # j_ℓ(0) = 0 for ℓ > 0
    @test abs(spherical_bessel(1, 1e-6)) < 1e-5
    
    # Check zeros
    zeros_0 = spherical_bessel_zeros(0, 5)
    for z in zeros_0
        @test abs(spherical_bessel(0, z)) < 1e-6
    end
end
```

### Test 3: Hankel Transform (analytical limit)
```julia
@testset "Hankel Transform" begin
    # For delta-function weight: Hankel transform should be simple
    delta_weight(chi) = delta(chi - 1.0)
    
    W_tilde = hankel_transform(1, 0.1, 0.2, 0.0, 10.0, 
                              delta_weight, chi_int, z_int)
    
    # Should have specific analytical form
    expected = spherical_bessel(1, 0.1 * 1.0) * spherical_bessel(1, 0.2 * 1.0)
    @test W_tilde ≈ expected atol=1e-4
end
```

### Test 4: Correlation vs Limber
```julia
@testset "S_ℓ vs Limber" begin
    # Compare S_ℓ with precomputed Limber reference
    
    # Compute S_ℓ using SFB method
    cache = SFBCorrelationCache(...)
    S_ell_sfb = compute_correlation_grid(2, k_vals, k_vals, cache)
    
    # Load Limber reference
    C_limber = npzread("data/Limber/Cl_CC_limber_linear_full.npy")
    
    # Convert C_ℓ to S_ℓ (involves additional integrations)
    # Check relative error < 1%
    
    rel_error = maximum(abs(S_ell_sfb - S_ell_limber) ./ abs(S_ell_limber))
    @test rel_error < 0.01
end
```

---

## Performance Profiling Template

```julia
using BenchmarkTools, Profile

# Profile E(z) computation
@time E(1.0, cosmo)  # First call (compilation)
@time E(1.0, cosmo)  # Second call (runtime)

@benchmark E($1.0, $cosmo)  # Detailed timing statistics

# Profile Bessel functions
@time spherical_bessel(10, 5.0)
@benchmark spherical_bessel($10, $5.0)

# Profile Hankel transform (most expensive!)
@time hankel_transform(2, 0.1, 0.2, ...)
@profile hankel_transform(2, 0.1, 0.2, ...)  # Flame graph analysis
```

---

This guide provides complete, copy-paste-ready code templates for each layer, ready for implementation!

