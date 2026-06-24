using Base.Threads
using LinearAlgebra

function get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number)
    z = FastTransforms.clenshawcurtisnodes(Float64, N)
    z = (zmax - zmin) / 2 * z .+ (zmin + zmax) / 2

    z[1] *= (1 - 1e-8)
    z[end] *= (1 + 1e-8) # TODO: replace with a principled endpoint treatment.

    return z
end

function get_clencurt_weights_z(zmin::Number, zmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (zmax - zmin) / 2 * w

    return w
end

# function bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
#                                  z_range::AbstractArray, n_cheb::Int, N::Number, chi_of_z::Any)
function bessel_cheb_eval_beyond(ℓ::Number, zmin::Number, zmax::Number, kmin::Number, kmax::Number, 
                                 z_range::AbstractArray, n_cheb::Int, N::Integer, chi_of_z::Any)
    Nz = length(z_range)
    k = get_clencurt_grid_z(kmin, kmax, N)
    z = get_clencurt_grid_z(zmin, zmax, N)

    # 1. Calcolo ottimizzato della matrice di Chebyshev T tramite ricorrenza
    T = zeros(n_cheb + 1, N)
    # Mappa z dall'intervallo [zmin, zmax] a [-1, 1]
    x = @. (2 * z - (zmax + zmin)) / (zmax - zmin)
    
    T[1, :] .= 1.0
    if n_cheb >= 1
        T[2, :] .= x
    end
    for i in 3:(n_cheb + 1)
        @. T[i, :] = 2 * x * T[i-1, :] - T[i-2, :]
    end

    # 2. Calcolo ottimizzato di Bessel (Zero allocazioni nel loop e Column-Major)
    Bessel = zeros(Nz, N)
    chi_vals = chi_of_z.(z_range) # Pre-calcolo per evitare di chiamare la funzione Nz * N volte

    # Parallelizziamo sul ciclo esterno (colonne) per sfruttare il column-major
    Threads.@threads for j in 1:N
        kj = k[j]
        for i in 1:Nz
            @inbounds Bessel[i, j] = sphericalbesselj(ℓ, chi_vals[i] * kj)
        end
    end

    return T, Bessel
end

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


