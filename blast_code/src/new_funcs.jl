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
function compute_T̃_beyond(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
    if kmin >= kmax 
        throw(DomainError("The integration range is unphysical. Make sure kmin < kmax.")) 
    end
    
    nχ = length(χ)
    nR = length(R)

    x = get_clencurt_grid(kmin, kmax, N)
    w = get_clencurt_weights(kmin, kmax, N)
    T, Bessel1 = bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)

    # Output: T_tilde[1, i, ridx, l] where:
    # i = index for χ
    # ridx = index for R (which parametrizes k' and k'')
    # l = Chebyshev polynomial index
    T_tilde = zeros(1, nχ, nR, n_cheb+1)
    
    for (ridx, r) in enumerate(R)
        # Bessel1[i,k] = jₗ(kχᵢ) - already computed
        # Bessel2[i,k] = jₗ(k'χᵢ) where k' is parametrized by R
        # Bessel3[i,k] = jₗ(k'χᵢ) - same as Bessel2 for the first χ₁ argument
        # Bessel4[i,k] = jₗ(k''χᵢ) where k'' is parametrized by R
        
        Bessel2 = zeros(nχ, N)  # jₗ(k'χᵢ)
        Bessel3 = zeros(nχ, N)  # jₗ(k'χᵢ) - for χ₁ in second window function
        Bessel4 = zeros(nχ, N)  # jₗ(k''χᵢ) - for χ₂ in second window function
        
        Threads.@threads for i in 1:nχ
            # For the Hankel-transformed window functions:
            # W̃ᴬ(k',k) contributes: jₗ(kχ₁) jₗ(k'χ₁)  → Bessel1 and Bessel2
            # W̃ᴮ(k'',k) contributes: jₗ(kχ₂) jₗ(k''χ₂) → Bessel1 (different χ) and Bessel4
            # But since R = χ₁/χ₂, we parametrize as:
            # k' → r*χ[i]  (for first window function)
            # k'' → r*χ[i]  (for second window function, but different χ argument)
            
            Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)  # jₗ(k'χ₁)
            Bessel3[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)  # jₗ(k'χ₁) for W̃ᴬ
            Bessel4[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)  # jₗ(k''χ₂)
        end

        α = w .* (x .^ β)  # Integration weight: f^AB(k) = k^β * w
        
        # Single integration over k with four Bessel function terms
        @tturbo l in 1:n_cheb+1, i in 1:nχ
                Cij = zero(eltype(w))
                for k in 1:N
                    # ∫ dk f(k) Tₙ(k) jₗ(kχ₁) jₗ(kχ₂) jₗ(k'χ₁) jₗ(k''χ₂)
                    # α[k] contains the k-dependent weight and integration measure
                    # T[l,k] is the Chebyshev polynomial (or power spectrum approximation)
                    # Bessel1[i,k], Bessel2[i,k], Bessel3[i,k], Bessel4[i,k] are the four Bessel terms
                    Cij += α[k] * T[l,k] * Bessel1[i,k] * Bessel2[i,k] * Bessel3[i,k] * Bessel4[i,k]
                end
                T_tilde[1,i,ridx,l] = Cij
    end

    return T_tilde

end