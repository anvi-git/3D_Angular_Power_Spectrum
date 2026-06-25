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

"""
    bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number,
                            z_range::AbstractArray, n_cheb::Int, N::Integer, chi_of_z::Any)
                            -> Tuple{Matrix{Float64}, Matrix{Float64}}
Build the Chebyshev basis matrix in `z` and the spherical-Bessel evaluation
matrix needed for a mixed Chebyshev/Bessel expansion.
The function:
1. constructs Clenshaw–Curtis grids in `k` and `z`,
2. evaluates Chebyshev polynomials \\(T\_n(x)\\) on the mapped `z` nodes,
3. evaluates `sphericalbesselj(ℓ, χ(z_i) k_j)` on the supplied redshift range.
# Arguments
- `ℓ::Number`: spherical Bessel order.
- `zmin::Number`: lower bound of the redshift interval used for the Chebyshev grid.
- `zmax::Number`: upper bound of the redshift interval used for the Chebyshev grid.
- `kmin::Number`: lower bound of the wavenumber interval used for the Clenshaw–Curtis grid.
- `kmax::Number`: upper bound of the wavenumber interval used for the Clenshaw–Curtis grid.
- `z_range::AbstractArray`: redshift samples where the Bessel kernel is evaluated.
- `n_cheb::Int`: highest Chebyshev order to include.
- `N::Integer`: number of Clenshaw–Curtis nodes.
- `chi_of_z::Any`: function mapping redshift `z` to comoving distance `χ(z)`.
# Returns
- `T::Matrix{Float64}`: Chebyshev basis matrix with size `(n_cheb + 1, N)`.
- `Bessel::Matrix{Float64}`: matrix with size `(length(z_range), N)` containing
  `sphericalbesselj(ℓ, chi_of_z(z_i) * k_j)`.
"""
function bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                                 z_range::AbstractArray, n_cheb::Int, N::Integer, chi_of_z::Any)
    Nz = length(z_range)
    k = get_clencurt_grid_z(kmin, kmax, N)
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

    Bessel = zeros(Nz, N)
    chi_vals = chi_of_z.(z_range)

    Threads.@threads for j in 1:N
        kj = k[j]
        for i in 1:Nz
            @inbounds Bessel[i, j] = sphericalbesselj(ℓ, chi_vals[i] * kj)
        end
    end

    return T, Bessel
end

# function nuova_bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
#                                  z_range::AbstractArray, n_cheb::Int, N::Integer, chi_of_z::Any)
#     Nz = length(z_range)
#     k = get_clencurt_grid_z(kmin, kmax, Nz)
#     z = get_clencurt_grid_z(zmin, zmax, N)

#     T = zeros(n_cheb + 1, N)
#     x = @. (2 * z - (zmax + zmin)) / (zmax - zmin)
#     T[1, :] .= 1.0
#     if n_cheb >= 1
#         T[2, :] .= x
#     end
#     for i in 3:(n_cheb + 1)
#         @. T[i, :] = 2 * x * T[i-1, :] - T[i-2, :]
#     end

#     Bessel = zeros(Nz, N)
#     chi_vals = chi_of_z.(z)

#     Threads.@threads for j in 1:Nz
#         kj = k[j]
#         for i in 1:N
#             @inbounds Bessel[j, i] = sphericalbesselj(ℓ, chi_vals[i] * kj)
#         end
#     end

#     return T, Bessel
# end

"""
    compute_W_tilde(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number,
                    z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)
                    -> Array{Float64,4}

Compute the 4D tensor `W_tilde` used in the Chebyshev/Bessel projection
pipeline.
The function constructs a Clenshaw–Curtis grid and weights, evaluates the
Chebyshev basis and spherical-Bessel kernels, and then accumulates the
weighted contraction explicitly over the quadrature index. The result is
stored in a tensor of shape `(1, Nz, Nz, n_cheb + 1)`.
# Arguments
- `ℓ::Number`: spherical Bessel order.
- `zmin::Number`: lower bound of the integration interval.
- `zmax::Number`: upper bound of the integration interval.
- `kmin::Number`: lower bound of the wavenumber interval.
- `kmax::Number`: upper bound of the wavenumber interval.
- `z_range::AbstractArray`: redshift samples used to evaluate `χ(z)`.
- `n_cheb::Int`: highest Chebyshev order to include.
- `N::Number`: number of Clenshaw–Curtis nodes.
- `chi_of_z::Any`: function mapping redshift `z` to comoving distance `χ(z)`.
# Returns
- `Array{Float64,4}`: tensor `W_tilde` with dimensions `(1, Nz, Nz, n_cheb + 1)`.
# Notes
- `chi = chi_of_z.(z_range)` is the comoving-distance sampling on the redshift grid.
- The current implementation assumes `Bessel2 == Bessel1`; if a second Bessel
  kernel is introduced, the contraction structure will need to be generalized.
- The current explicit assignment fills the second and third tensor dimensions
  with the same contracted value for each `p`, so this layout is meant as an
  intermediate storage form and may be optimized later.
"""
# function compute_W_tilde(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
#                          z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)
#     if zmin >= zmax 
#         throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
#     end

#     Nz = length(z_range)
#     chi = chi_of_z.(z_range)
#     k = get_clencurt_grid_z(kmin, kmax, N)
#     w = get_clencurt_weights_z(zmin, zmax, N)    
#     T, Bessel1 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
#     T_tilde = zeros(1, Nz, Nz, n_cheb+1)
    
#     #when Bessel2 = Bessel 1: comment these four lines below
#     # Bessel2 = zeros(Nz, N)        
#     # Threads.@threads for i in 1:Nz
#     #     Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
#     # end
#     α = w 

#     #commenting the 6 lines of code below works fine with N = 2^5 +1, with N=2^15 is 5 steps in 3 minutes
#     # for (p, chi_val) in enumerate(chi) 
#     #     Bessel2 = zeros(Nz, N)        
#     #     Threads.@threads for i in 1:Nz
#     #         Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
#     #     end
#     #     α = w #β = 2 for CC, -2 for LL and 0 for CL.
#     for (p, chi_val) in enumerate(chi)
#         @tturbo for l in 1:n_cheb+1, i in 1:Nz
#             Cij = zero(eltype(w))
#             for k in 1:N
#                 #Cij +=  T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
#                 Cij +=  T[l,k] * Bessel1[i,k] * Bessel1[i,k] * α[k] #when Bessel2 = Bessel1
#             end
#             T_tilde[1,i,p,l] = Cij
#         end
#     end

#     return T_tilde

# end

"""
    compute_W̃(ℓ::Number, zmin::Real, zmax::Real, kmin::Real, kmax::Real,
              z_range::AbstractArray, n_cheb::Int, N::Int, chi_of_z::Any)
              -> Array{Float64,4}
Compute the intermediate kernel tensor `W̃` used in the Chebyshev/Bessel
projection pipeline.
The function first builds the Clenshaw–Curtis quadrature weights on
`[zmin, zmax]`, evaluates the Chebyshev basis and spherical-Bessel matrix via
`bessel_cheb_eval_beyond`, and then combines them into a precomputed weighted
matrix. The final result is expanded into a 4D tensor with shape
`(1, Nk, Nk, n_cheb + 1)`.
- `A` is formed by weighting the squared Bessel kernel column-wise with the
  quadrature weights.
- The current implementation fills each `T_tilde[1, i, p, l]` with the same
  value `C[i, l]` for all `p`, so the third dimension is currently redundant.
# Arguments
- `ℓ::Number`: spherical Bessel order.
- `zmin::Real`: lower bound of the integration interval.
- `zmax::Real`: upper bound of the integration interval.
- `kmin::Real`: lower bound of the wavenumber interval.
- `kmax::Real`: upper bound of the wavenumber interval.
- `z_range::AbstractArray`: redshift samples used for the Bessel evaluation.
- `n_cheb::Int`: highest Chebyshev order to include.
- `N::Int`: number of Clenshaw–Curtis nodes.
- `chi_of_z::Any`: function mapping redshift `z` to comoving distance `χ(z)`.
# Returns
- `Array{Float64,4}`: a tensor `T_tilde` with dimensions `(1, Nk, Nk, n_cheb + 1)`.
"""
function compute_W̃(ℓ::Number, zmin::Real, zmax::Real, kmin::Real, kmax::Real, 
                         z_range::AbstractArray, n_cheb::Int, N::Int, chi_of_z::Any)
    if zmin >= zmax 
        throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
    end

    Nk = length(z_range)
    chi = chi_of_z.(z_range)
    w = get_clencurt_weights_z(zmin, zmax, N)    
    T, Bessel1 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
    
    # if Bessel2 is different from Bessel1, then
    # Bessel2 = zeros(Nk, N)
    #_, Bessel2 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
    #A = @. Bessel1 * Bessel2 * w'
    # Pre-computation of the weight matrix: (Nk x N)
    # Multiply every column of Bessel1.^2 for the corresponding weight w[k]
    # w' transforms the vector into a row matrix (1 x N) for correct broadcasting
    A = @. Bessel1^2 * w' 
    # A is (Nk x N), T' is (N x n_cheb+1) -> C will be (Nk x n_cheb+1)
    C = A * T'
    T_tilde = zeros(1, Nk, Nk, n_cheb+1)    
    for l in 1:n_cheb+1
        for p in 1:Nk
            for i in 1:Nk
                @inbounds T_tilde[1, i, p, l] = C[i, l]
            end
        end
    end

    return T_tilde
end


