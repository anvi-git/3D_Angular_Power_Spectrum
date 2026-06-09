using Base.Threads
using LinearAlgebra

########## common to the two versions of compute_W_tilde ##########
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

####### version 1 #######
function bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                                 z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)

    Nz = length(z_range)
    k = get_clencurt_grid_z(kmin, kmax, N)
    z = get_clencurt_grid_z(zmin, zmax, N)
    z_cheb = chebpoints(n_cheb, zmin,zmax) 
    c = FastChebInterp.ChebPoly(z_cheb, SA[zmin], SA[zmax])
    chi = chi_of_z.(z_cheb)

    T = zeros(n_cheb+1,N) 
    Threads.@threads for i in 1:n_cheb+1
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.
        T[i,:] = copy_c.(z)
    end

    Bessel = zeros(Nz, N)
    Threads.@threads for i in 1:Nz
            #Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
            Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi_of_z(z_range[i]) * k)
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
####### end of version 1 #######
############################
####### version 2 #######
function build_chebyshev_matrix(z::AbstractVector, zmin::Number, zmax::Number, n_cheb::Int)
    N = length(z)
    z_norm = @. 2 * (z - zmin) / (zmax - zmin) - 1
    z_norm = clamp.(z_norm, -1.0, 1.0)

    T = Matrix{Float64}(undef, N, n_cheb + 1)
    T[:, 1] .= 1.0
    if n_cheb >= 1
        T[:, 2] .= z_norm
    end
    for l in 3:n_cheb+1
        @. T[:, l] = 2 * z_norm * T[:, l-1] - T[:, l-2]
    end
    return T  # (N, n_cheb+1)
end

function precompute_bessel_matrices(ℓ::Number, chi::AbstractVector, k::AbstractVector)
    Nz = length(chi)
    N  = length(k)
    B = Matrix{Float64}(undef, Nz, N)

    Threads.@threads for i in 1:Nz
        @views B[i, :] = SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k)
    end
    return B  # (Nz, N)
end

function compute_W_tilde_opt(ℓ::Number, zmin::Number, zmax::Number,
                                            kmin::Number, kmax::Number,
                                            z_range::AbstractArray, n_cheb::Int,
                                            N::Number, chi_of_z::Any)

    if zmin >= zmax
        throw(DomainError("The integration range is unphysical. Make sure zmin < zmax."))
    end

    Nz = length(z_range)

    k = get_clencurt_grid_z(kmin, kmax, N)
    z = get_clencurt_grid_z(zmin, zmax, N)

    # stessi pesi del file 1
    α = get_clencurt_weights_z(zmin, zmax, N)

    # stessa T del file 1, ma costruita in layout BLAS-friendly
    T = build_chebyshev_matrix(z, zmin, zmax, n_cheb)  # (N, n_cheb+1)

    # nel file 1 Bessel1 viene da chi_of_z.(z_cheb), ma la dimensione allocata dipende da Nz;
    # per mantenere la semantica effettiva usata in compute_W_tilde, qui usiamo chi(z_range)
    chi = chi_of_z.(z_range)
    Bessel1 = precompute_bessel_matrices(ℓ, chi, k)  # (Nz, N)

    # nel file 1 Bessel2 è identica a Bessel1 ad ogni iterazione su p
    BW = similar(Bessel1)
    @. BW = Bessel1 * Bessel1 * α'   # broadcast riga per riga con α

    # contrazione su k: (Nz,N) * (N,n_cheb+1) -> (Nz,n_cheb+1)
    slice_ipℓ = Matrix{Float64}(undef, Nz, n_cheb + 1)
    mul!(slice_ipℓ, BW, T)

    # stessa shape del file 1
    T_tilde = Array{Float64}(undef, 1, Nz, Nz, n_cheb + 1)

    # nel file 1 ogni p produce lo stesso slice
    for p in 1:Nz
        @views T_tilde[1, :, p, :] .= slice_ipℓ
    end

    return T_tilde
end
####### end of version 2 #######