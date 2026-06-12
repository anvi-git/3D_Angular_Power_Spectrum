using Base.Threads
using LinearAlgebra

"""
    get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number)

Construct a Clenshaw–Curtis grid in redshift over the interval `[zmin, zmax]`.

The function:
1. generates `N` Clenshaw–Curtis nodes on `[-1, 1]`,
2. linearly maps them to `[zmin, zmax]`,
3. applies a tiny perturbation to the first and last nodes to avoid exact endpoint
   coincidences that can cause downstream numerical issues in integrations or
   interpolations.

# Arguments
- `zmin::Number`: Lower redshift bound.
- `zmax::Number`: Upper redshift bound.
- `N::Number`: Number of quadrature/grid nodes.

# Returns
- `z::Vector{Float64}`: Clenshaw–Curtis nodes mapped to the redshift interval.

# Notes
- The endpoint adjustment is a temporary workaround; a more principled treatment
  should be preferred if exact endpoints matter for the application.
"""
function get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number)
    z = FastTransforms.clenshawcurtisnodes(Float64, N)
    z = (zmax - zmin) / 2 * z .+ (zmin + zmax) / 2

    z[1] *= (1 - 1e-8)
    z[end] *= (1 + 1e-8) # TODO: replace with a principled endpoint treatment.

    return z
end

"""
    get_clencurt_weights_z(zmin::Number, zmax::Number, N::Number)

Construct Clenshaw–Curtis quadrature weights over the redshift interval `[zmin, zmax]`.

The function first builds the Clenshaw–Curtis quadrature object for `N` points,
computes the corresponding weights on the reference interval `[-1, 1]`, and then
rescales them to the physical interval `[zmin, zmax]`.

# Arguments
- `zmin::Number`: Lower redshift bound.
- `zmax::Number`: Upper redshift bound.
- `N::Number`: Number of quadrature points.

# Returns
- `w::Vector{Float64}`: Clenshaw–Curtis quadrature weights rescaled to `[zmin, zmax]`.

# Notes
- The factor `(zmax - zmin) / 2` comes from the linear map between `[-1, 1]` and
  `[zmin, zmax]`.
- These weights are intended to be used together with the nodes returned by
  `get_clencurt_grid_z`.
"""
function get_clencurt_weights_z(zmin::Number, zmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (zmax - zmin) / 2 * w

    return w
end

"""
    bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number,
                            kmin::Number, kmax::Number,
                            z_range::AbstractArray, n_cheb::Int,
                            N::Number, chi_of_z::Any)

Precompute Chebyshev basis functions and spherical Bessel kernels for the
beyond-Limber integration.

The function constructs Clenshaw–Curtis grids in `k` and `z`, builds the
Chebyshev basis over the interval `[zmin, zmax]`, evaluates each basis function
on the `z` integration grid, and computes the spherical Bessel functions
`j_ℓ(k χ(z))` on the redshift sampling given by `z_range`.

# Arguments
- `ℓ::Number`: Multipole order of the spherical Bessel function.
- `zmin::Number`: Lower bound of the redshift interval used for the Chebyshev expansion.
- `zmax::Number`: Upper bound of the redshift interval used for the Chebyshev expansion.
- `kmin::Number`: Minimum wavenumber of the Clenshaw–Curtis integration grid.
- `kmax::Number`: Maximum wavenumber of the Clenshaw–Curtis integration grid.
- `z_range::AbstractArray`: Redshift values where the Bessel kernel is evaluated.
- `n_cheb::Int`: Order of the Chebyshev expansion.
- `N::Number`: Number of Clenshaw–Curtis quadrature points.
- `chi_of_z::Any`: Function or callable object returning the comoving distance `χ(z)`.

# Returns
- `T::Matrix{Float64}`: Matrix of size `(n_cheb + 1, N)` containing the Chebyshev
  basis functions evaluated on the Clenshaw–Curtis `z` grid.
- `Bessel::Matrix{Float64}`: Matrix of size `(length(z_range), N)` containing
  `j_ℓ(k χ(z))` evaluated on the Clenshaw–Curtis `k` grid for each value in `z_range`.

# Notes
- The rows of `T` correspond to individual Chebyshev basis polynomials.
- The rows of `Bessel` correspond to different redshift samples in `z_range`.
- The computation is threaded both for the Chebyshev basis construction and for
  the Bessel kernel evaluation.
- The variable `chi = chi_of_z.(z_cheb)` is currently computed but not used in the
  final Bessel evaluation, so it may be removable unless needed for future changes.
"""
function bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                                 z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)

    Nz = length(z_range)
    k = get_clencurt_grid_z(kmin, kmax, N)
    z = get_clencurt_grid_z(zmin, zmax, N)
    z_cheb = chebpoints(n_cheb, zmin, zmax) 
    c = FastChebInterp.ChebPoly(z_cheb, SA[zmin], SA[zmax])
    #chi = chi_of_z.(z_cheb)

    T = zeros(n_cheb + 1, N) 
    Threads.@threads for i in 1:n_cheb+1
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.
        T[i, :] = copy_c.(z)
    end

    Bessel = zeros(Nz, N)
    Threads.@threads for i in 1:Nz
        #Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
        Bessel[i, :] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi_of_z(z_range[i]) * k)
    end

    return T, Bessel
end

"""
    compute_W_tilde(ℓ::Number, zmin::Number, zmax::Number,
                    kmin::Number, kmax::Number,
                    z_range::AbstractArray, n_cheb::Int,
                    N::Number, chi_of_z::Any)

Compute the pre-integrated kernel `W̃` used in the beyond-Limber projected-matter
density calculation.

This function first checks that the redshift interval is physical, then builds the
Clenshaw–Curtis grids and weights, evaluates the Chebyshev basis and spherical
Bessel kernels, and finally contracts them to form the precomputed tensor
`T_tilde`. In the current implementation, the Bessel kernel is evaluated with
`Bessel1` reused in place of `Bessel2`, so the routine effectively computes the
auto-kernel case `j_ℓ(kχ)^2` weighted by the Clenshaw–Curtis quadrature rule.

# Arguments
- `ℓ::Number`: Multipole order of the spherical Bessel function.
- `zmin::Number`: Lower bound of the redshift interval.
- `zmax::Number`: Upper bound of the redshift interval.
- `kmin::Number`: Minimum wavenumber for the Clenshaw–Curtis `k` grid.
- `kmax::Number`: Maximum wavenumber for the Clenshaw–Curtis `k` grid.
- `z_range::AbstractArray`: Redshift samples at which the kernel is evaluated.
- `n_cheb::Int`: Number of Chebyshev polynomials used in the expansion.
- `N::Number`: Number of Clenshaw–Curtis nodes.
- `chi_of_z::Any`: Callable mapping redshift `z` to comoving distance `χ(z)`.

# Returns
- `T_tilde::Array{Float64,4}`: A 4D tensor of shape `(1, Nz, Nz, n_cheb + 1)`
  containing the quadrature-contracted Chebyshev/Bessel kernel.

# Notes
- The input range must satisfy `zmin < zmax`; otherwise a `DomainError` is thrown.
- `χ(z)` is precomputed on `z_range` for the outer loop and the same Bessel block
  is used on both sides of the contraction in the current auto-correlation setup.
- The commented `Bessel2` block indicates where the general two-kernel case would
  enter if `Bessel1` and `Bessel2` were different.
- The tensor is designed for the later construction of the projected matter density
  `w̃` in the non-Limber angular power spectrum pipeline.
"""
function compute_W_tilde(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                         z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)
    if zmin >= zmax 
        throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
    end

    Nz = length(z_range)
    chi = chi_of_z.(z_range)
    k = get_clencurt_grid_z(kmin, kmax, N)
    w = get_clencurt_weights_z(zmin, zmax, N)    
    T, Bessel1 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
    T_tilde = zeros(1, Nz, Nz, n_cheb+1)
    
# when Bessel2 = Bessel 1 I comment these four lines below
    # Bessel2 = zeros(Nz, N)        
    # Threads.@threads for i in 1:Nz
    #     Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
    # end
    α = w #β = 2 for CC, -2 for LL and 0 for CL.

#commenting the 6 lines of code below works fine with N = 2^5 +1, with N=2^15 is 5 steps in 3 minutes
    # for (p, chi_val) in enumerate(chi) 
    #     Bessel2 = zeros(Nz, N)        
    #     Threads.@threads for i in 1:Nz
    #         Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
    #     end
    #     α = w #β = 2 for CC, -2 for LL and 0 for CL.
    for (p, chi_val) in enumerate(chi)
        @tturbo for l in 1:n_cheb+1, i in 1:Nz
            Cij = zero(eltype(w))
            for k in 1:N
                #Cij +=  T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
                Cij +=  T[l,k] * Bessel1[i,k] * Bessel1[i,k] * α[k] #when Bessel2 = Bessel1
            end
            T_tilde[1,i,p,l] = Cij
        end
    end

    return T_tilde

end