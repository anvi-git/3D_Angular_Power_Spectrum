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
    compute_T̃_general(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number, k1::Number, k2::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
Compute the generalized spherical kernel by reusing `compute_T̃` and dressing the result with the extra external Bessel factors.

# Arguments
- `ℓ::Number`: Multipole order.
- `χ::AbstractArray`: Array containing values of the comoving distance.
- `R::AbstractArray`: Array containing values for the `R=χ₂/χ₁` variable.
- `kmin::Number` and `kmax::Number`: Integration range in k.
- `β::Number`: Exponent of the k dependence of the integral.
- `k1::Number`: External wavenumber associated with the first leg.
- `k2::Number`: External wavenumber associated with the second leg.
- `n_cheb::Int`: Number of Chebyshev polynomials used in the approximation of the power spectra.
- `N::Int`: Number of integration points in k.
"""
function compute_T̃_general(
    ℓ::Number,
    χ::AbstractArray,
    R::AbstractArray,
    kmin::Number,
    kmax::Number,
    β::Number,
    k1::Number,
    k2::Number;
    n_cheb::Int = 119,
    N::Int = 2^(15)+1,
)
    if kmin >= kmax
        throw(DomainError("The integration range is unphysical. Make sure kmin < kmax."))
    end

    nχ = length(χ)
    nR = length(R)

    x = get_clencurt_grid(kmin, kmax, N)
    w = get_clencurt_weights(kmin, kmax, N)
    T, Bessel1 = bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)

    J1 = SpecialFunctions.sphericalbesselj.(ℓ, k1 .* χ)
    J2 = zeros(nχ, nR)
    for (ridx, r) in enumerate(R)
        J2[:, ridx] = SpecialFunctions.sphericalbesselj.(ℓ, k2 .* χ .* r)
    end

    T_tilde = zeros(1, nχ, nR, n_cheb + 1)

    for (ridx, _) in enumerate(R)
        Bessel2 = zeros(nχ, N)

        Threads.@threads for i in 1:nχ
            Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, R[ridx] * χ[i] * x)
        end

        α = w .* (x .^ β) #β = 2 for CC, -2 for LL and 0 for CL.

        for l in 1:n_cheb+1, i in 1:nχ
            Cij = zero(eltype(w))
            for k in 1:N
                Cij += T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
            end
            T_tilde[1,i,ridx,l] = J1[i] * J2[i,ridx] * Cij
        end
    end

    return T_tilde
end