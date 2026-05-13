"""
    window_prefactor(Hz, nz; bias=one.(Hz), growth=one.(Hz), c=Blast.C_LIGHT)

Tracer-dependent prefactor f^A(χ) used in the generalized window function.
"""
window_prefactor(Hz, nz; bias=one.(Hz), growth=one.(Hz), c=Blast.C_LIGHT) = (Hz ./ c) .* bias .* nz .* growth

"""
    generalized_window_function(k, χ; ell::Integer, weight, r = identity)

Generalized window W_l^A(k,χ) = f^A(χ) j_l(k r(χ)).
"""
function generalized_window_function(k, χ; ell::Integer, weight, r = identity)
    radial = r === identity ? χ : r.(χ)
    return weight .* SpecialFunctions.sphericalbesselj.(ell, k .* radial)
end

"""
    compute_T̃_beyond(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
Compute integrals of the Bessel functions and the Chebyshev polynomials with four Bessel function arguments.
This is the precomputation part of the code for the generalized BLAST formalism with Hankel-transformed window functions.

Computes:
    T̃ₙ;ℓᴬᴮ(χ₁, χ₂) ≡ ∫_{kmin}^{kmax} dk f^{AB}(k) T_n(k) j_ℓ(kχ₁) j_ℓ(kχ₂) j_ℓ(k'χ₁) j_ℓ(k''χ₂)

where k', k'' are wavenumber arguments of the Hankel-transformed window functions W̃ᴬ(k',k) and W̃ᴮ(k'',k).

# Arguments
- `ℓ::Number`: Multipole order

- `χ::AbstractArray`: Array containing values of the comoving distance. 

- `R::AbstractArray`: Array containing values for the R=χ₁/χ₂ variable.

- `kmin::Number` and `kmax::Number`: Integration range in k.

- `β::Number`: Exponent of the k dependence of the integral. This parameter depends on the combination of tracers: β=2,-2,0 for clustering, cosmic shear and the cross-correlation respectively.

- `n_cheb::Int`: Number of chebyshev polynomials used in the approximation of the power spectra.

- `N::Int`: Number of integration points in k.
"""
function compute_T̃_beyond(ℓ::Number, χ::AbstractArray, R::AbstractArray, 
                           kp::AbstractArray, kpp::AbstractArray,  # k' and k'' grids
                           kmin::Number, kmax::Number, β::Number; 
                           n_cheb::Int = 119, N::Int = 2^(15)+1)
    
    nχ = length(χ)
    nR = length(R)
    nkp = length(kp)
    nkpp = length(kpp)

    x = Blast.get_clencurt_grid(kmin, kmax, N)
    w = Blast.get_clencurt_weights(kmin, kmax, N)
    T, Bessel_k = Blast.bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)

    # kp and kpp are passed as arguments, already in linear space
    # (not reassigned here - would cause log-space multiplication error)
    # Output: (nχ, nR, nkp, nkpp, n_cheb+1)
    T_tilde = zeros(nχ, nR, nkp, nkpp, n_cheb+1)
    
    α = w .* (x .^ β)
    
    for i in 1:nχ, ridx in 1:nR
        χ₁ = χ[i]
        χ₂ = χ₁ / R[ridx]  # or χ₁ * R[ridx], depending on convention
        
        # Bessel factors for k: j_ℓ(k·χ₁) and j_ℓ(k·χ₂)
        Bessel_k_χ1 = @views Bessel_k[i, :]  # Pre-computed j_ℓ(k·χ₁)
        Bessel_k_χ2 = zeros(N)
        Threads.@threads for n in 1:N
            Bessel_k_χ2[n] = sphericalbesselj(ℓ, x[n] * χ₂)  # j_ℓ(k·χ₂)
        end
        
        # Loop over fixed k' and k''
        for kpi in 1:nkp, kppi in 1:nkpp
            # These DON'T depend on k (fixed arguments)
            bessel_kp_χ1 = sphericalbesselj(ℓ, kp[kpi] * χ₁)
            bessel_kpp_χ2 = sphericalbesselj(ℓ, kpp[kppi] * χ₂)
            
            # Integrate over k with Chebyshev coefficients
            for l in 1:n_cheb+1
                integral = 0.0
                @simd for n in 1:N
                    integral += α[n] * T[l, n] * 
                               Bessel_k_χ1[n] * Bessel_k_χ2[n] * bessel_kp_χ1 * bessel_kpp_χ2
                end
                T_tilde[i, ridx, kpi, kppi, l] = integral
            end
        end
    end

    return T_tilde
end
