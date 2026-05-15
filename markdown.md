# Hankel Transform Extension for BLAST Code: Implementation Plan

**Project Goal**: Extend BLAST.jl (arXiv:2410.03632) to compute 3D power spectra with Hankel-transformed window functions using spherical Bessel decomposition. Result: 4 Bessel function arguments (2 existing + 2 new from Hankel transforms).

---

## Core Equation

$$S_\ell^{\text{den}}(k_1,k_2) = N_\ell^{\text{den}} \int dk \, k^2 P(k) \, \tilde{W}_\ell^{\text{den}}(k_1,k) \, \tilde{W}_\ell^{\text{den}}(k_2,k)$$

where the Hankel-transformed window is:

$$\tilde{W}_\ell^{\text{den}}(k_1,k) = \int d\chi \, \chi^2 \, f^{\text{den}}(\chi) \, j_\ell(k\chi) \, j_\ell(k_1\chi)$$

**Four Bessel functions**:
1. $j_\ell(k \cdot \chi)$ — radial Bessel in power spectrum integral
2. $j_\ell(k_1 \cdot \chi)$ — Hankel transform first argument
3. $j_\ell(k \cdot \chi)$ — appears again in second window
4. $j_\ell(k_2 \cdot \chi)$ — Hankel transform second argument

---

## Implementation Steps (In Order)

### STEP 1: Spherical Bessel Utilities
**File**: Create `src/spherical_bessel_utils.jl`

```julia
function spherical_bessel_zeros(ℓ::Int, n_zeros::Int)
    # Return first n_zeros zeros of j_ℓ(x)
    # Use SpecialFunctions.sphericalbesselj with root-finding
end

function spherical_bessel_derivative(ℓ::Int, x::AbstractArray)
    # Compute d/dx[j_ℓ(x)] using recurrence: j'_ℓ = j_{ℓ-1} - (ℓ+1)/x * j_ℓ
    return sphericalbesselj.(ℓ-1, x) .- (ℓ+1)./x .* sphericalbesselj.(ℓ, x)
end

function spherical_bessel_grid(ℓ::Int, χ_min::Number, χ_max::Number, n_grid::Int)
    # Create adaptive grid denser where j_ℓ oscillates
end
```

### STEP 2: Window Prefactor
**File**: Create `src/window_prefactor.jl`

```julia
struct WindowPrefactor
    H_z::Vector{Float64}      # Hubble parameter H(z)
    bias_z::Vector{Float64}   # Galaxy bias b(z)
    nz::Vector{Float64}       # Number density n(z)
    growth_z::Vector{Float64} # Growth factor D(z)
    χ::Vector{Float64}        # Comoving distances
end

function compute_window_prefactor(z::AbstractArray, cosmo, bias_z, nz, growth_z)
    H_z = Blast.compute_hubble_factor.(z, Ref(cosmo))
    χ = Blast.compute_χ.(z, Ref(cosmo))
    return WindowPrefactor(H_z, bias_z, nz, growth_z, χ)
end
```

### STEP 3: Hankel Window Kernel
**File**: Create `src/hankel_window_kernel.jl`

```julia
function hankel_window_kernel(ℓ::Int, k_out::Number, k::Number, 
                               χ::AbstractArray, prefactor::WindowPrefactor; 
                               rtol=1e-6, c=Blast.C_LIGHT)
    # Computes: ∫ dχ χ² f^den(χ) j_ℓ(kχ) j_ℓ(k_out*χ)
    
    function integrand(χ_val)
        # Interpolate prefactor quantities at χ_val
        f_den = (H(χ_val)/c) * bias(χ_val) * n(χ_val) * growth(χ_val)
        return χ_val^2 * f_den * sphericalbesselj(ℓ, k*χ_val) * 
               sphericalbesselj(ℓ, k_out*χ_val)
    end
    
    result, _ = quadgk(integrand, first(χ), last(χ), rtol=rtol)
    return result
end

function compute_hankel_window_grid(ℓ::Int, k_grid::AbstractArray, 
                                     k_out_grid::AbstractArray, 
                                     χ::AbstractArray, 
                                     prefactor::WindowPrefactor; rtol=1e-6)
    # Output: W_tilde[i_kout, i_k] shape (length(k_out_grid), length(k_grid))
    
    n_kout = length(k_out_grid)
    n_k = length(k_grid)
    W_tilde = zeros(n_kout, n_k)
    
    Threads.@threads for i_kout in 1:n_kout
        k_out = k_out_grid[i_kout]
        for i_k in 1:n_k
            k = k_grid[i_k]
            W_tilde[i_kout, i_k] = hankel_window_kernel(ℓ, k_out, k, χ, prefactor; rtol=rtol)
        end
    end
    
    return W_tilde
end

function compute_hankel_window_pair(ℓ::Int, k1_grid::AbstractArray, 
                                     k2_grid::AbstractArray, χ::AbstractArray, 
                                     prefactor_A::WindowPrefactor, 
                                     prefactor_B::WindowPrefactor; rtol=1e-6)
    # Auto-correlation: compute once and alias
    if prefactor_A === prefactor_B
        W = compute_hankel_window_grid(ℓ, k1_grid, k2_grid, χ, prefactor_A; rtol=rtol)
        return W, W
    else
        # Cross-correlation: compute both
        W_A = compute_hankel_window_grid(ℓ, k1_grid, k2_grid, χ, prefactor_A; rtol=rtol)
        W_B = compute_hankel_window_grid(ℓ, k1_grid, k2_grid, χ, prefactor_B; rtol=rtol)
        return W_A, W_B
    end
end

function compute_T̃_hankel(ℓ::Int, χ::AbstractArray, k_grid::AbstractArray, 
                           k_out_grid::AbstractArray, kmin::Number, kmax::Number, 
                           prefactor::WindowPrefactor; 
                           n_cheb::Int=119, N::Int=2^15+1)
    # Output: T_hankel[i_kout, i_χ, i_cheb] 
    # shape: (length(k_out_grid), length(χ), n_cheb+1)
    
    nχ = length(χ)
    nk_out = length(k_out_grid)
    
    # 1. Clenshaw-Curtis grid and weights
    x = Blast.get_clencurt_grid(kmin, kmax, N)
    w = Blast.get_clencurt_weights(kmin, kmax, N)
    
    # 2. Chebyshev polynomials + first set of Bessel functions
    T, Bessel_k = Blast.bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)
    
    # 3. Hankel-transformed window functions
    W_hankel = compute_hankel_window_grid(ℓ, x, k_out_grid, χ, prefactor)
    # Shape: (length(k_out_grid), length(x))
    
    # 4. Main contraction
    T_hankel = zeros(nk_out, nχ, n_cheb + 1)
    
    @tturbo for l in 1:n_cheb+1, i_χ in 1:nχ, i_k_out in 1:nk_out
        Cij = zero(eltype(w))
        for n in 1:N
            Cij += T[l, n] * Bessel_k[i_χ, n] * W_hankel[i_k_out, n] * w[n]
        end
        T_hankel[i_k_out, i_χ, l] = Cij
    end
    
    return T_hankel
end

function compute_Sl_hankel(ℓ::Int, k_out_grid::AbstractArray, χ::AbstractArray, 
                            T_hankel::AbstractArray, cheb_coeff::AbstractArray; 
                            normalize=true)
    # Output: S_ℓ[i_k1, i_k2] shape (length(k_out_grid), length(k_out_grid))
    
    S_ℓ = @tullio S[i1, i2] := T_hankel[i1, j, l] * T_hankel[i2, j, l] * cheb_coeff[j, l]
    
    if normalize
        S_ℓ *= (2/π)^2  # SFB normalization
    end
    
    return S_ℓ
end

include("spherical_bessel_utils.jl")

export hankel_window_kernel, compute_hankel_window_grid, compute_hankel_window_pair
export compute_T̃_hankel, compute_Sl_hankel, WindowPrefactor, compute_window_prefactor
export spherical_bessel_zeros, spherical_bessel_derivative, spherical_bessel_grid

using Test, Blast

@testset "Hankel Transform Suite" begin
    
    @testset "Hankel Kernel" begin
        W_tilde = compute_hankel_window_grid(2, k_grid, k_out_grid, χ, prefactor)
        @test size(W_tilde) == (length(k_out_grid), length(k_grid))
        @test all(isfinite.(W_tilde))
        @test norm(W_tilde) > 0
    end
    
    @testset "T̃ with Hankel" begin
        T_hankel = compute_T̃_hankel(2, χ, k_grid, k_out_grid, kmin, kmax, prefactor)
        @test size(T_hankel) == (length(k_out_grid), length(χ), n_cheb+1)
        @test all(isfinite.(T_hankel))
    end
    
    @testset "Power Spectrum" begin
        S_ℓ = compute_Sl_hankel(2, k_out_grid, χ, T_hankel, cheb_coeff)
        @test size(S_ℓ) == (length(k_out_grid), length(k_out_grid))
        @test all(isfinite.(S_ℓ))
    end
    
    @testset "Limber Limit" begin
        # For large ℓ, 4-Bessel should approach Limber (2-Bessel)
        # Implement comparison test
    end
end

# 1. Setup
ℓ = 2  # or vector
k_grid = exp.(range(log(0.01), log(10), length=100))
k_out_grid = exp.(range(log(0.01), log(10), length=50))
χ = LinRange(26, 7000, 30)

# 2. Build prefactor (from existing background data)
z_array = z_of_χ.(χ)
prefactor = Blast.compute_window_prefactor(z_array, cosmo, bias_z, growth_z)

# 3. Compute Hankel transforms (expensive; do once per ℓ)
W_tilde = Blast.compute_hankel_window_grid(ℓ, k_grid, k_out_grid, χ, prefactor)

# 4. Compute T̃ coefficients
T_hankel = Blast.compute_T̃_hankel(ℓ, χ, k_grid, k_out_grid, kmin, kmax, prefactor, 
                                    n_cheb=119, N=2^15+1)

# 5. Chebyshev coefficients of power spectrum (existing code)
cheb_coeff = zeros(length(χ), n_cheb + 1)
vals = power_spectrum.(10 .^ k_cheb, χ[1], χ[1])
plan = Blast.plan_fft(vals)
for j in 1:length(χ)
    vals = power_spectrum.(10 .^ k_cheb, χ[j], χ[j])
    cheb_coeff[j, :] = Blast.fast_chebcoefs(vals, plan)
end

# 6. Final contraction
S_ℓ = Blast.compute_Sl_hankel(ℓ, k_out_grid, χ, T_hankel, cheb_coeff)

# 7. Validate against Limber
# Compare with existing BLAST results