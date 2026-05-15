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


"""
    computehankelwindowCC(ℓ::Int, χ::AbstractArray, k_window::AbstractArray,
                          W_g::AbstractArray) → WtildeCC::Array{Float64,2}

Compute the Hankel-transformed galaxy-clustering (CC) window function

    W̃_ℓ^{den}(k, χ₀) = ∫ dχ' χ'² W_g(χ') j_ℓ(k χ₀) j_ℓ(k χ')

which is the auto-correlation (A = B = den) specialisation of `computehankelwindow`.

# Arguments
- `ℓ::Int`                  Multipole order
- `χ::AbstractArray`        Comoving distance grid (length n)
- `k_window::AbstractArray` Wavenumber grid for the Hankel convolution (length n_k)
- `W_g::AbstractArray`      Real-space galaxy window W_g(χ) evaluated on χ (length n)

# Returns
- `WtildeCC::Array{Float64,2}` shape (n_k, n)
"""
function computehankelwindowCC(ℓ::Int, χ::AbstractArray,
                                k_window::AbstractArray,
                                W_g::AbstractArray)
    return computehankelwindow(ℓ, χ, k_window, W_g)
end


"""
    computewCCwithwindows(ℓ::Int, χ::AbstractArray, R::AbstractArray,
                          kmin::Number, kmax::Number,
                          Wtilde1::AbstractArray, Wtilde2::AbstractArray,
                          k_window::AbstractArray,
                          cheb_coeff::AbstractArray;
                          n_cheb::Int = 119, N::Int = 2151)
        → w_CC::Matrix{Float64}

Full pipeline: given pre-computed Hankel windows `Wtilde1` and `Wtilde2`
(each of shape n_k × n), build the beyond-BLAST coefficient tensor T̃_CC and
immediately contract it with the Chebyshev-decomposed power spectrum to
produce w_CC(χ, R) for a single multipole ℓ.

# Arguments
- `ℓ::Int`                   Multipole order
- `χ::AbstractArray`         Comoving distance grid (length n)
- `R::AbstractArray`         Chebyshev R-grid (length nR)
- `kmin,kmax::Number`        Integration range in k
- `Wtilde1,Wtilde2`          Hankel-transformed windows, shape (n_k, n)
- `k_window::AbstractArray`  Wavenumber grid matching rows of Wtilde (length n_k)
- `cheb_coeff::AbstractArray` Chebyshev coefficients of P(k,χ,R), shape (n, nR, n_cheb+1)
- `n_cheb::Int`              Number of Chebyshev polynomials
- `N::Int`                   Integration points in k

# Returns
- `w_CC::Matrix{Float64}` shape (n, nR)
"""
function computewCCwithwindows(ℓ::Int, χ::AbstractArray, R::AbstractArray,
                                kmin::Number, kmax::Number,
                                Wtilde1::AbstractArray, Wtilde2::AbstractArray,
                                k_window::AbstractArray,
                                cheb_coeff::AbstractArray;
                                n_cheb::Int = 119, N::Int = 2151)
    # Build T̃_CC  →  shape (1, n, nR, n_cheb+1)
    Ttilde = computeTwithwindows(ℓ, χ, R, kmin, kmax,
                                  Wtilde1, Wtilde2, k_window;
                                  n_cheb = n_cheb, N = N)
    # Drop the leading ℓ-dimension  →  (n, nR, n_cheb+1)
    Ttilde_single = dropdims(Ttilde; dims = 1)
    # Chebyshev contraction  →  (n, nR)
    w_CC = dropdims(sum(cheb_coeff .* Ttilde_single; dims = 3); dims = 3)
    return w_CC
end


"""
    compute_T̃_with_windows(ℓ::Number, χ::AbstractArray, R::AbstractArray, 
                           kmin::Number, kmax::Number, 
                           W_tilde_1::AbstractArray, W_tilde_2::AbstractArray,
                           k_window::AbstractArray;
                           n_cheb::Int = 119, N::Int = 2^(15)+1)

Compute BLAST coefficients with Hankel-transformed window functions.

This is the "proper beyond-BLAST" for galaxy clustering with realistic windows.

# Arguments
- `ℓ::Number`: Multipole order
- `χ::AbstractArray`: Comoving distance array
- `R::AbstractArray`: R = χ₁/χ₂ array
- `kmin, kmax::Number`: Integration range in k
- `W_tilde_1, W_tilde_2::AbstractArray`: Hankel-transformed windows [nk × nχ]
- `k_window::AbstractArray`: Wavenumber grid for windows

# Returns
- `T_tilde::Array{Float64, 4}`: Precomputed coefficients [1, nχ, nR, n_cheb+1]
"""
function compute_T̃_with_windows(ℓ::Number, χ::AbstractArray, R::AbstractArray, 
                                kmin::Number, kmax::Number, 
                                W_tilde_1::AbstractArray, W_tilde_2::AbstractArray,
                                k_window::AbstractArray;
                                n_cheb::Int = 119, N::Int = 2^(15)+1)
    
    nχ = length(χ)
    nR = length(R)
    nk_win = length(k_window)
    
    # Get Clenshaw-Curtis grid and weights
    x = get_clencurt_grid(kmin, kmax, N)
    w = get_clencurt_weights(kmin, kmax, N)
    
    # Compute Chebyshev polynomials and Bessel functions
    T, Bessel1 = bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)
    
    T_tilde = zeros(1, nχ, nR, n_cheb+1)
    
    for (ridx, r) in enumerate(R)
        Bessel2 = zeros(nχ, N)
        
        # Compute second set of Bessel functions
        Threads.@threads for i in 1:nχ
            Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)
        end
        
        # Weight array (β = 2 for clustering)
        α = w .* (x .^ 2)
        
        # Loop over window wavenumbers
        for i_k1 in 1:nk_win
            k1 = k_window[i_k1]
            
            for i_k2 in 1:nk_win
                k2 = k_window[i_k2]
                
                # Additional Bessel functions from windows
                for i in 1:nχ
                    i_χ2 = Int(round(i / r))  # Index for χ₂ = χ₁/r
                    if i_χ2 < 1 || i_χ2 > nχ
                        continue
                    end
                    
                    # Window contributions
                    W1_factor = W_tilde_1[i_k1, i]
                    W2_factor = W_tilde_2[i_k2, i_χ2]
                    
                    # Main integration loop
                    for l in 1:n_cheb+1
                        Cij = 0.0
                        for k in 1:N
                            # Standard BLAST part
                            base = T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
                            
                            # Window modification
                            Cij += base * W1_factor * W2_factor
                        end
                        T_tilde[1, i, ridx, l] += Cij
                    end
                end
            end
        end
    end
    
    return T_tilde
end