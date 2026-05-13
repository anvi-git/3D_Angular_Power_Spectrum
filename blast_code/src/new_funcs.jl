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
                               kp::AbstractArray, kpp::AbstractArray, 
                               kmin::Number, kmax::Number;  
                               n_cheb::Int = 119, N::Int = 2^(15)+1)
    nχ = length(χ)
    nR = length(R)
    nkp = length(kp)
    nkpp = length(kpp)

    x   = Blast.get_clencurt_grid(kmin, kmax, N)
    w   = Blast.get_clencurt_weights(kmin, kmax, N)
    T, Bessel_k = Blast.bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)

    # Precompute the three weight arrays once, outside all loops
    α_m2 = w .* (x .^ (-2))   # β = -2  (LL)
    α_0  = w                   # β =  0  (CL), x^0 = 1
    α_p2 = w .* (x .^ 2)      # β = +2  (CC)

    T_LL = zeros(nχ, nR, nkp, nkpp, n_cheb+1)
    T_CL = zeros(nχ, nR, nkp, nkpp, n_cheb+1)
    T_CC = zeros(nχ, nR, nkp, nkpp, n_cheb+1)

    Threads.@threads for i in 1:nχ
        Bessel_k_χ2 = zeros(N)

        for ridx in 1:nR
            χ₁ = χ[i]
            χ₂ = χ₁ / R[ridx]
            Bessel_k_χ1 = @views Bessel_k[i, :]

            for n in 1:N
                Bessel_k_χ2[n] = sphericalbesselj(ℓ, x[n] * χ₂)
            end

            for kpi in 1:nkp, kppi in 1:nkpp
                bessel_kp_χ1  = sphericalbesselj(ℓ, kp[kpi]  * χ₁)
                bessel_kpp_χ2 = sphericalbesselj(ℓ, kpp[kppi] * χ₂)
                const_factor  = bessel_kp_χ1 * bessel_kpp_χ2

                for l in 1:n_cheb+1
                    I_m2 = 0.0; I_0 = 0.0; I_p2 = 0.0
                    @simd for n in 1:N
                        base = T[l, n] * Bessel_k_χ1[n] * Bessel_k_χ2[n]
                        I_m2 += α_m2[n] * base
                        I_0  += α_0[n]  * base
                        I_p2 += α_p2[n] * base
                    end
                    T_LL[i, ridx, kpi, kppi, l] = I_m2 * const_factor
                    T_CL[i, ridx, kpi, kppi, l] = I_0  * const_factor
                    T_CC[i, ridx, kpi, kppi, l] = I_p2 * const_factor
                end
            end
        end
    end

    return T_LL, T_CL, T_CC
end