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

function W_tilde_computed(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number,
                             Nk::Int, Nkp::Int, n_cheb::Int, N::Number, chi_of_z::Any)

    if zmin >= zmax
        throw(DomainError("The integration range is unphysical. Make sure zmin < zmax."))
    end

    kp = get_clencurt_grid_z(kmin, kmax, Nkp)
    w  = get_clencurt_weights_z(zmin, zmax, N)

    T, Bessel1 = bessel_func(ℓ, zmin, zmax, kmin, kmax, Nk, n_cheb, N, chi_of_z)
    T_tilde = zeros(eltype(w), Nk, Nkp, n_cheb + 1, 1)
    _, Bessel2 = bessel_func(ℓ, zmin, zmax, kmin, kmax, Nkp, n_cheb, N, chi_of_z)

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

function bessel_func(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                                 Nk::Int, n_cheb::Int, N::Integer, chi_of_z::Any)
    
    k = get_clencurt_grid_z(kmin, kmax, Nk)
    z = get_clencurt_grid_z(zmin, zmax, N)

    T = zeros(n_cheb + 1, N)
    x = @. (2 * z - (zmax + zmin)) / (zmax - zmin)
    T[1, :] .= 1.0
    if n_cheb >= 1
        T[2, :] .= x
    end
    for i in 3:(n_cheb + 1)
        @. T[i, :] = 2 * x * T[i-1, :] - T[i-2, :]
    end
    
    Bessel = zeros(Nk, N)
    chi_vals = chi_of_z.(z)

    Threads.@threads for j in 1:Nk
        kj = k[j]
        for i in 1:N
#           Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, χhi[i] * x)
            Bessel[j, i] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi_vals[i] * kj)
            #@inbounds Bessel[j, i] = sphericalbesselj(ℓ, chi_vals[i] * kj)
        end
    end

    return T, Bessel

end


# function compute_W̃(ℓ::Number, zmin::Real, zmax::Real, kmin::Real, kmax::Real, 
#                          z_range::AbstractArray, n_cheb::Int, N::Int, chi_of_z::Any)
#     if zmin >= zmax 
#         throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
#     end

#     Nk = length(z_range)
#     chi = chi_of_z.(z_range)
#     w = get_clencurt_weights_z(zmin, zmax, N)    
#     T, Bessel1 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
    
#     # if Bessel2 is different from Bessel1, then
#     # Bessel2 = zeros(Nk, N)
#     #_, Bessel2 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
#     #A = @. Bessel1 * Bessel2 * w'
#     # Pre-computation of the weight matrix: (Nk x N)
#     # Multiply every column of Bessel1.^2 for the corresponding weight w[k]
#     # w' transforms the vector into a row matrix (1 x N) for correct broadcasting
# #    A = @. Bessel1^2 * w' #this squares Bessel1 and multiplies it by the conjugate transpose of w
#     A = @. Bessel1 * Bessel1 * w' #this squares Bessel1 and multiplies it by the conjugate transpose of w
#     # A is (Nk x N), T' is (N x n_cheb+1) -> C will be (Nk x n_cheb+1)
#     C = A * T'
#     T_tilde = zeros(1, Nk, Nk, n_cheb+1)    
#     for l in 1:n_cheb+1
#         for p in 1:Nk
#             for i in 1:Nk
#                 @inbounds T_tilde[1, i, p, l] = C[i, l]
#             end
#         end
#     end

#     return T_tilde
# end