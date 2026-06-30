using Base.Threads
using LinearAlgebra

"""
    get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number) -> Vector{Float64}
Construct a Clenshaw–Curtis grid in the redshift interval `[zmin, zmax]`.
The nodes are first generated on the interval `[-1, 1]` using
`FastTransforms.clenshawcurtisnodes`, then mapped to `[zmin, zmax]`.
# Arguments
- `zmin::Number`: lower bound of the redshift interval.
- `zmax::Number`: upper bound of the redshift interval.
- `N::Number`: number of Clenshaw–Curtis nodes to generate.
# Returns
- `Vector{Float64}`: the mapped redshift grid.
"""
function get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number)
    z = FastTransforms.clenshawcurtisnodes(Float64, N)
    z = (zmax - zmin) / 2 * z .+ (zmin + zmax) / 2

    z[1] *= (1 - 1e-8)
    z[end] *= (1 + 1e-8) # TODO: replace with a principled endpoint treatment.

    return z
end

"""
    get_clencurt_weights_z(zmin::Number, zmax::Number, N::Number) -> Vector{Float64}
Construct Clenshaw–Curtis quadrature weights mapped to the redshift interval `[zmin, zmax]`.
The weights are first computed on the canonical Clenshaw–Curtis grid using
`FastTransforms.chebyshevmoments1` together with `FastTransforms.clenshawcurtisweights`, 
then rescaled to the physical interval by the factor `(zmax - zmin) / 2`.
# Arguments
- `zmin::Number`: lower bound of the redshift interval.
- `zmax::Number`: upper bound of the redshift interval.
- `N::Number`: number of quadrature weights.
# Returns
- `Vector{Float64}`: quadrature weights on `[zmin, zmax]`.
"""
function get_clencurt_weights_z(zmin::Number, zmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (zmax - zmin) / 2 * w

    return w
end

function W_tilde_computation(ℓ::Number, x_min::Number, x_max::Number, kmin::Number, kmax::Number,
                             Nk::Int, Nkp::Int, n_cheb::Int, N::Number, x::AbstractVector)

    if x_min >= x_max
        throw(DomainError("The integration range is unphysical. Make sure x_min < x_max."))
    end

    kp = get_clencurt_grid_z(kmin, kmax, Nkp)
    w  = get_clencurt_weights_z(x_min, x_max, N)

    T, Bessel1 = bessel_computation(ℓ, x_min, x_max, kmin, kmax, Nk, n_cheb, N, x)
    T_tilde = zeros(eltype(w), Nk, Nkp, n_cheb + 1, 1)
    _, Bessel2 = bessel_computation(ℓ, x_min, x_max, kmin, kmax, Nkp, n_cheb, N, x)

    for ic in 1:(n_cheb + 1) 
        @tturbo for ik in 1:Nk, ikp in 1:Nkp
            Cij = zero(eltype(w))
            for iN in 1:N
                Cij += T[ic, iN] * Bessel1[ik, iN] * Bessel2[ikp, iN] * w[iN] * kp[ikp]
            end
        T_tilde[ik, ikp, ic, 1] = Cij
        end
    end
    return T_tilde

end

function bessel_computation(ℓ::Number, x_min::Number, x_max::Number, kmin::Number, kmax::Number, 
                                 Nk::Int, n_cheb::Int, N::Integer, x::AbstractVector)
    
    k = get_clencurt_grid_z(kmin, kmax, Nk)
    x_grid = get_clencurt_grid_z(x_min, x_max, N)

    T = zeros(n_cheb + 1, N)
    xx = @. (2 * x_grid - (x_max + x_min)) / (x_max - x_min)
    T[1, :] .= 1.0
    if n_cheb >= 1
        T[2, :] .= xx
    end
    for i in 3:(n_cheb + 1)
        @. T[i, :] = 2 * xx * T[i-1, :] - T[i-2, :]
    end
    
    Bessel = zeros(Nk, N)

    Threads.@threads for j in 1:Nk
        kj = k[j]
        for i in 1:N
            #here x[i] is the i-th element of the x array
            Bessel[j, i] = @views SpecialFunctions.sphericalbesselj.(ℓ, x[i] * kj) 
        end
    end

    return T, Bessel

end