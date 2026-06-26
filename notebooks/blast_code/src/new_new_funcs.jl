function compute_W_again(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
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
    α = w 
    for (p, chi_val) in enumerate(chi)
        @tturbo for l in 1:n_cheb+1, i in 1:Nz
            Cij = zero(eltype(w))
            for ik in 1:N
                Cij +=  T[l,ik] * Bessel1[i,ik] * Bessel1[i,ik] * α[ik] #when Bessel2 = Bessel1
            end
            T_tilde[1,i,p,l] = Cij
        end
    end

    return T_tilde

end