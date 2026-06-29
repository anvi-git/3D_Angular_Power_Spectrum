# function new_compute_W_tilde(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
#                          z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)
#     if zmin >= zmax 
#         throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
#     end

#     Nz = length(z_range)
#     Nk = Nz
#     Nkp = Nz - 20
#     kp = get_clencurt_grid_z(kmin, kmax, Nkp)
#     k = get_clencurt_grid_z(kmin, kmax, Nk)
#     z = get_clencurt_grid_z(zmin, zmax, N)
#     w = get_clencurt_weights_z(zmin, zmax, N)    
#     T, Bessel1 = bessel_cheb_eval_beyond(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)
#     T_tilde = zeros(1, Nk, Nkp, n_cheb+1)
#     Bessel2 = zeros(N, Nkp) 

#     chi_vals = chi_of_z.(z)
#     Threads.@threads for j in 1:Nkp
#         kj = kp[j]
#         for i in 1:N
#             @inbounds Bessel2[i,j] = sphericalbesselj(ℓ, chi_vals[i] * kj)
#     end

#     α = w

#     @tturbo for ic in 1:n_cheb+1, ik in 1:Nk, ikp in 1:Nkp
#         Cij = zero(eltype(w))
#         for iN in 1:N
#             Cij +=  T[ic,iN] * Bessel1[iN,ik] * Bessel1[iN,ikp] * α[iN] #when Bessel2 = Bessel1
#         end
#         T_tilde[1,ik,ikp,ic] = Cij
#     end

#     return T_tilde

# end

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