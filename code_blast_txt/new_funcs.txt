using Base.Threads

function get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    z = FastTransforms.clenshawcurtisnodes(Float64, N)
    z = (zmax - zmin) / 2 * z .+ (zmin + zmax) / 2 

    z[1] *= (1-1e-8)
    z[end] *= (1+1e-8) #TODO: this is just a quick patch, need to figure this out properly.

    return z
end

function get_clencurt_weights_z(zmin::Number, zmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (zmax - zmin) / 2 * w

    return w
end

function bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                                 z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)

    Nz = length(z_range)
    k = get_clencurt_grid_z(kmin, kmax, N)
    z = get_clencurt_grid_z(zmin, zmax, N)
    z_cheb = chebpoints(n_cheb, zmin,zmax) 
    c = FastChebInterp.ChebPoly(z_cheb, SA[zmin], SA[zmax])
    chi = chi_of_z.(z_cheb) # calcolare chi(z) su z_range

    T = zeros(n_cheb+1,N) 
    Threads.@threads for i in 1:n_cheb+1
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.
        T[i,:] = copy_c.(z)
    end

    Bessel = zeros(Nz, N)
    Threads.@threads for i in 1:Nz
            Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
    end

    return T, Bessel

end

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
    
    for (p, chi_val) in enumerate(chi)
        Bessel2 = zeros(Nz, N)        
        Threads.@threads for i in 1:Nz
            Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
        end
        α = w #β = 2 for CC, -2 for LL and 0 for CL.
         
        @tturbo for l in 1:n_cheb+1, i in 1:Nz
            Cij = zero(eltype(w))
            for k in 1:N
                Cij +=  T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
            end
            T_tilde[1,i,p,l] = Cij
        end
    end

    return T_tilde

end

