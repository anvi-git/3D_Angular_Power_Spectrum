function get_clencurt_grid_z(zmin::Number, zmax::Number, N::Number)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    x = (zmax - zmin) / 2 * x .+ (zmin + zmax) / 2 
    
    x[1] *= (1 - 1e-8)
    x[end] *= (1 + 1e-8) # Patch dal tuo file originale [cite: 2]
    return x
end

function generate_cheb_matrix_z(z_nodes, zmin::Number, zmax::Number)
    nz = length(z_nodes)
    
    # Inizializziamo il polinomio di Chebyshev nell'intervallo z lineare
    # (Nota: qui usiamo la scala lineare, non log10 come per k, perché z varia di poco)
    c = FastChebInterp.ChebPoly(z_nodes, SA[zmin], SA[zmax])
    
    # Matrice di trasformazione (ordine_polinomio x punti_griglia)
    T_z = zeros(nz, nz) 
    
    Threads.@threads for i in 1:nz
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.0
        # Valutiamo il polinomio i-esimo su tutti i nodi di z
        T_z[i, :] = copy_c.(z_nodes)
    end
    return T_z
end

function compute_geometric_tensor(ℓ::Number, zmin::Number, zmax::Number, z_cheb_nodes::AbstractVector, x_array::AbstractVector, k_grid::AbstractVector, k_ext::AbstractVector, n_z::Integer; N::Integer = N)
    Nk = length(k_grid)
    Next = length(k_ext)
    
    # 1. Griglia fine di integrazione in z e relativi pesi di Clenshaw-Curtis
    z_int_nodes = get_clencurt_grid_z(zmin, zmax, N)
    z_int_weights = FastTransforms.clenshawcurtisweights(FastTransforms.chebyshevmoments1(Float64, N)) * (zmax - zmin) / 2
    
    interp_χ = LinearInterpolation(z_cheb_nodes, x_array, extrapolation = ExtrapolationType.Linear)
    χ_int = interp_χ.(z_int_nodes)
    
    T_z = generate_cheb_matrix_z(z_int_nodes, zmin, zmax)
    
    # NOVITÀ: Includiamo i pesi di Clenshaw-Curtis direttamente nei polinomi
    T_weighted = T_z .* z_int_weights'
    
    # 3. PRE-CALCOLO DELLE BESSEL SU k_grid (Risparmia milioni di valutazioni)
    # Creiamo una matrice [N, Nk]
    Bessel_v = zeros(N, Nk)
    Threads.@threads for i in 1:Nk
        k_v = k_grid[i]
        for n in 1:N
            Bessel_v[n, i] = sphericalbesselj(ℓ, k_v * χ_int[n])
        end
    end
    
    # 4. Allocazione del Tensore Geometrico G [n_z, Nk, Next]
    G = zeros(n_z, Nk, Next)
    
    # Buffer pre-allocati per ogni thread per avere ZERO allocazioni nel ciclo
    n_threads = Threads.nthreads()
    bessel_e_buffers = [zeros(N) for _ in 1:n_threads]
    M_buffers = [zeros(n_z, N) for _ in 1:n_threads]
    
    # 5. Ciclo parallelo su Next
    # 5. Build Bessel_e as a true [N, Next] matrix
    Bessel_e = zeros(N, Next)
    Threads.@threads for j in 1:Next
        k_e = k_ext[j]
        @inbounds for n in 1:N
            Bessel_e[n, j] = sphericalbesselj(ℓ, k_e * χ_int[n])
        end
    end

    @tullio G[m, i, j] := T_weighted[m, n] * Bessel_e[n, j] * Bessel_v[n, i]

    return G
end
# function compute_geometric_tensor(ℓ::Number, zmin::Number, zmax::Number, z_cheb_nodes::AbstractVector, x_array::AbstractVector, k_grid::AbstractVector, k_ext::AbstractVector, n_z::Integer; N::Integer = N)
#     Nk = length(k_grid)
#     Next = length(k_ext)
    
#     # 1. Griglia fine di integrazione in z e relativi pesi di Clenshaw-Curtis
#     z_int_nodes = get_clencurt_grid_z(zmin, zmax, N)
#     z_int_weights = FastTransforms.clenshawcurtisweights(FastTransforms.chebyshevmoments1(Float64, N)) * (zmax - zmin) / 2
    
#     interp_χ = LinearInterpolation(z_cheb_nodes, x_array, extrapolation = ExtrapolationType.Linear)
#     χ_int = interp_χ.(z_int_nodes)
#     # 2. Pre-calcoliamo i polinomi di Chebyshev della base valutati sulla griglia fine di integrazione
#     # Usiamo lo stesso principio del tuo file originale
#     c_base = FastChebInterp.ChebPoly(get_clencurt_grid_z(zmin, zmax, n_z), SA[zmin], SA[zmax])
#     T_polys = zeros(n_z, N)
#     for m in 1:n_z
#         copy_c = deepcopy(c_base)
#         copy_c.coefs .*= 0
#         copy_c.coefs[m] = 1.0
#         T_polys[m, :] = copy_c.(z_int_nodes)
#     end
    
#     # 3. Allocazione del Tensore Geometrico G [n_z, Nk, Next]
#     # G[m, i, j] conterrà l'integrale del polinomio m-esimo con Bessel(k_i) e Bessel(k_ext_j)
#     G = zeros(n_z, Nk, Next)
    
#     # Ciclo parallelo per riempire il tensore geometrico
#     Threads.@threads for j in 1:Next
#         k_e = k_ext[j]
#         for i in 1:Nk
#             k_v = k_grid[i]
            
#             # Calcoliamo l'integrand geometrico combinando le due Bessel e i pesi di CC
#             # j_ℓ(k * χ) * j_ℓ(k_ext * χ) * dz_weights
#             bessel_prod_weights = @. sphericalbesselj(ℓ, k_v * χ_int) * sphericalbesselj(ℓ, k_e * χ_int) * z_int_weights
            
#             for m in 1:n_z
#                 # Contrazione lungo la griglia di integrazione z (prodotto scalare)
#                 G[m, i, j] = sum(T_polys[m, :] .* bessel_prod_weights)
#             end
#         end
#     end
    
#     return G
# end
######################################################
######################################################


# Definisci una volta sola fuori dai loop
# struct IntegrationGrid
#     z::Vector{Float64}
#     w::Vector{Float64}
#     jac::Float64
# end

# # Inizializza una volta sola
# function setup_grid(zmin, zmax, N)
#     CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
#     t_cc   = FastTransforms.clenshawcurtisnodes(Float64, N)
#     w_cc   = FastTransforms.clenshawcurtisweights(CC_obj)
#     z_cc   = @. (zmax + zmin) / 2 + (zmax - zmin) / 2 * t_cc
#     return IntegrationGrid(z_cc, w_cc, (zmax - zmin) / 2)
# end

# function bessel_cheb_eval(ℓ::Number, chi::AbstractArray, zed::AbstractArray, kmin::Number, kmax::Number, zmin::Number, zmax::Number; n_cheb::Int, N::Number)

#     chi_z = DataInterpolations.AkimaInterpolation(chi, zed)
#     n_chi = length(chi)
#     k = get_clencurt_grid(kmin, kmax, N)
#     z = get_clencurt_grid(zmin, zmax, N)
#     z_cheb = chebpoints(n_cheb, log10(zmin), log10(zmax)) 
#     c = FastChebInterp.ChebPoly(z_cheb, SA[log10(zmin)], SA[log10(zmax)])

#     W = zeros(n_cheb+1,N) 
#     Threads.@threads for i in 1:n_cheb+1
#         copy_c = deepcopy(c) 
#         copy_c.coefs .*= 0 
#         copy_c.coefs[i] = 1.
#         W[i,:] = copy_c.(log10.(z))
#     end

#     Bessel = zeros(n_chi, N)
#     Threads.@threads for i in 1:n_chi
#             Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi_z(zed[i]) * k)
#     end

#     return W, Bessel

# end

# ############################
# function bessel_cheb_eval(
#     grid::IntegrationGrid,
#     ℓ::Number,
#     chi::AbstractArray,
#     zed::AbstractArray,
#     kmin::Number,
#     kmax::Number;
#     n_cheb::Int
# )
#     chi_z = DataInterpolations.AkimaInterpolation(chi, zed)
#     n_chi = length(chi)
#     k = get_clencurt_grid(kmin, kmax, length(grid.w))
#     z_cheb = chebpoints(n_cheb, log10(first(grid.z)), log10(last(grid.z)))
#     c = FastChebInterp.ChebPoly(z_cheb, SA[log10(first(grid.z))], SA[log10(last(grid.z))])

#     W = zeros(n_cheb + 1, length(grid.z))
#     for i in 1:n_cheb+1
#         copy_c = deepcopy(c)
#         copy_c.coefs .= 0
#         copy_c.coefs[i] = 1.0
#         W[i, :] = copy_c.(log10.(grid.z))
#     end

#     Bessel = zeros(n_chi, length(k))
#     Threads.@threads for i in 1:n_chi
#         Bessel[i, :] = SpecialFunctions.sphericalbesselj.(ℓ, chi_z(zed[i]) * k)
#     end

#     return W, Bessel
# end

# function compute_W̃(ℓ::Number, chi::AbstractArray, zed::AbstractArray, kmin::Number, kmax::Number, zmin::Number, zmax::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
#     if zmin >= zmax 
#         throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
#     end

#     nz = length(zed)
#     n_chi = length(chi)
#     chi_z = DataInterpolations.AkimaInterpolation(chi, zed)
#     z = get_clencurt_grid(zmin, zmax, N)
#     w = get_clencurt_weights(kmin, kmax, N)
#     k = get_clencurt_grid(kmin, kmax, N)
#     kp = k
#     W, Bessel1 = bessel_cheb_eval(ℓ, chi, zed, kmin, kmax, zmin, zmax; n_cheb, N)

#     W_tilde = zeros(1, nz, nz, n_cheb+1)
    
#     for (p, q) in enumerate(chi)
#         Bessel2 = zeros(n_chi, N)
        
#         Threads.@threads for i in 1:n_chi
#             Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi_z(zed[i]) * kp)
#         end

#         α = w
         
#         @tturbo for l in 1:n_cheb+1, i in 1:nz
#             Cij = zero(eltype(w))
#             for k in 1:N
#                 Cij +=  W[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
#             end
#             W_tilde[1,i,p,l] = Cij
#         end
#     end

#     return W_tilde

# end

# function compute_W̃(
#     grid::IntegrationGrid,
#     ℓ::Number,
#     chi::AbstractRange,
#     zed::AbstractArray,
#     kmin::Number,
#     kmax::Number;
#     n_cheb::Int = 119
# )
#     if first(grid.z) <= last(grid.z)
#         throw(DomainError("The integration range is unphysical. Make sure zmin < zmax."))
#     end

#     n_z = length(zed)
#     n_chi = length(chi)
#     chi_z = DataInterpolations.AkimaInterpolation(chi, zed)
#     k = get_clencurt_grid(kmin, kmax, length(grid.w))
#     kp = k
#     W, Bessel1 = bessel_cheb_eval(grid, ℓ, chi, zed, kmin, kmax; n_cheb)
#     W_tilde = zeros(1, n_z, n_z, n_cheb+1)
    
#     for (p, q) in enumerate(zed)
#         Bessel2 = zeros(n_chi, length(grid.w))
        
#         Threads.@threads for i in 1:n_chi
#             Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, chi_z(zed[i]) * kp)
#         end

#         α = grid.w
         
#         @tturbo for l in 1:n_cheb+1, i in 1:n_z
#             Cij = zero(eltype(w))
#             for k in 1:length(grid.w)
#                 Cij +=  W[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
#             end
#             W_tilde[1,i,p,l] = Cij
#         end
#     end

#     return W_tilde

# end

# function compute_W̃_chebcoeffs(
#     chi_interp,       # interpolatore χ(z)
#     n_z,              # funzione n(z)
#     b_z,              # funzione b(z) (bias)
#     D_z,              # funzione D(z) (growth factor)
#     #H_z,              # funzione H(z)
#     zmin::Number,
#     zmax::Number;
#     n_cheb::Int,
#     N::Int = 2^15 + 1
# )
#     # Stessi nodi e pesi usati in compute_W_tilde_modes
#     CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
#     t_cc   = FastTransforms.clenshawcurtisnodes(Float64, N)
#     w_cc   = FastTransforms.clenshawcurtisweights(CC_obj)

#     z_cc = @. (zmax + zmin) / 2 + (zmax - zmin) / 2 * t_cc
#     χ_cc = chi_interp.(z_cc)
#     jac  = (zmax - zmin) / 2

#     # Valuta W_tilde su tutti i nodi CC
#     #W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc) * H_z(z_cc)
#     W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc)

#     # Costruisci T_mat: stessa logica di compute_W_tilde_modes
#     c_basis = FastChebInterp.ChebPoly(zeros(n_cheb+1), SA[-1.0], SA[1.0])
#     T_mat = zeros(n_cheb+1, N)
#     Threads.@threads for n in 1:n_cheb+1
#         basis = deepcopy(c_basis)
#         basis.coefs[n] = 1.0
#         T_mat[n, :] = basis.(t_cc)
#     end

#     # Proiezione: c_n = jac * sum_t w_cc[t] * W_vals[t] * T_n(t_cc[t])
#     c_n = zeros(n_cheb+1)
#     @tturbo for n in 1:n_cheb+1
#         s = zero(Float64)
#         for t in 1:N
#             s += w_cc[t] * W_vals[t] * T_mat[n, t]
#         end
#         c_n[n] = jac * s
#     end

#     return c_n  # vettore [n_cheb]
# end
###################
# function compute_W̃_chebcoeffs(
#     grid::IntegrationGrid,
#     chi_interp,
#     n_z,
#     b_z,
#     D_z;
#     n_cheb::Int = 119
# )
#     N = length(grid.z)
#     z_cc = grid.z
#     w_cc = grid.w
#     jac  = grid.jac

#     center = (first(z_cc) + last(z_cc)) / 2
#     t_cc = @. (z_cc - center) / jac

#     χ_cc = chi_interp.(z_cc)

#     W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc)

#     c_basis = FastChebInterp.ChebPoly(zeros(n_cheb+1), SA[-1.0], SA[1.0])
#     T_mat = zeros(n_cheb+1, N)
#     Threads.@threads for n in 1:n_cheb+1
#         basis = deepcopy(c_basis)
#         basis.coefs[n] = 1.0
#         T_mat[n, :] = basis.(t_cc)
#     end

#     c_n = zeros(n_cheb+1)
#     @tturbo for n in 1:n_cheb+1
#         s = zero(Float64)
#         for t in 1:N
#             s += w_cc[t] * W_vals[t] * T_mat[n, t]
#         end
#         c_n[n] = jac * s
#     end

#     return c_n
# end

# function compute_W̃_chebcoeffs(
#     chi_interp,       # interpolatore χ(z)
#     n_z,              # funzione n(z)
#     b_z,              # funzione b(z) (bias)
#     D_z,              # funzione D(z) (growth factor)
#     #H_z,              # funzione H(z)
#     zmin::Number,
#     zmax::Number;
#     n_cheb::Int,
#     N::Int = 2^15 + 1
# )
#     # Stessi nodi e pesi usati in compute_W_tilde_modes
#     CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
#     t_cc   = FastTransforms.clenshawcurtisnodes(Float64, N)
#     w_cc   = FastTransforms.clenshawcurtisweights(CC_obj)

#     z_cc = @. (zmax + zmin) / 2 + (zmax - zmin) / 2 * t_cc
#     χ_cc = chi_interp.(z_cc)
#     jac  = (zmax - zmin) / 2

#     # Valuta W_tilde su tutti i nodi CC
#     #W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc) * H_z(z_cc)
#     W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc)

#     # Costruisci T_mat: stessa logica di compute_W_tilde_modes
#     c_basis = FastChebInterp.ChebPoly(zeros(n_cheb+1), SA[-1.0], SA[1.0])
#     T_mat = zeros(n_cheb+1, N)
#     Threads.@threads for n in 1:n_cheb+1
#         basis = deepcopy(c_basis)
#         basis.coefs[n] = 1.0
#         T_mat[n, :] = basis.(t_cc)
#     end# # Proiezione: c_n = jac * sum_t w_cc[t] * W_vals[t] * T_n(t_cc[t])
#     c_n = zeros(n_cheb+1)
#     @tturbo for n in 1:n_cheb+1
#         s = zero(Float64)
#         for t in 1:N
#             s += w_cc[t] * W_vals[t] * T_mat[n, t]
#         end
#         c_n[n] = jac * s
#     end

#     return c_n  # vettore [n_cheb]
# end

# function compute_Sℓ(
#     I_A::AbstractArray,
#     I_B::AbstractArray,
#     k_array::AbstractVector;
#     Pk::AbstractVector = ones(length(k_array)),
#     normalization::Real = 2/pi
# )
#     nk = length(k_array)

#     w_k = Blast.simpson_weight_array(nk)
#     Δk = nk > 1 ? (last(k_array) - first(k_array)) / (nk - 1) : one(eltype(k_array))

#     integrand_k = normalization .* k_array.^2 .* Pk .* w_k .* Δk

#     @tullio S[li, ik1, ik2] := integrand_k[m] *
#         I_A[li, ik1, m] *
#         I_B[li, ik2, m]

#     return S
# end