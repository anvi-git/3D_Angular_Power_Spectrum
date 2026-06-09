using Base.Threads
using LinearAlgebra

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

####################################
##################################
#################################
# # ─────────────────────────────────────────────────────────────────────────────
# # ─────────────────────────────────────────────────────────────────────────────
# #  Costruzione matrice T dei polinomi di Chebyshev
# #
# #  OPT 3: eliminato deepcopy in loop threaded.
# #  T[k, i] = T_i(z_k), layout (N, n_cheb+1) colonna-major:
# #  il loop interno su k (riga) è ora cache-friendly.
# #
# #  OPT 4 (layout): T è ora (N, n_cheb+1) invece di (n_cheb+1, N),
# #  così T[k, l] nel loop BLAS è accesso colonna — nessun miss di cache.
# # ─────────────────────────────────────────────────────────────────────────────

function build_chebyshev_matrix(z::AbstractVector, zmin::Number, zmax::Number, n_cheb::Int)
    N = length(z)

    # Mappa i nodi z in [-1, 1] per la formula dei polinomi di Chebyshev
    z_norm = @. 2 * (z - zmin) / (zmax - zmin) - 1
    z_norm  = clamp.(z_norm, -1.0, 1.0)

    # T[k, l] = cos((l-1) * acos(z_norm[k])), con l = 1..n_cheb+1
    # Costruzione via ricorrenza a tre termini: stabile e O(N * n_cheb)
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

# # ─────────────────────────────────────────────────────────────────────────────
# #  Precomputa le due matrici di Bessel su griglie k potenzialmente diverse
# #
# #  OPT 1 (principale): Bessel1 e Bessel2 sono entrambe calcolate UNA SOLA
# #  VOLTA qui, fuori da qualsiasi loop su p. Costo: O(Nz * N) invece di
# #  O(Nz² * N). Con Nz=100 questo elimina 99% delle valutazioni di Bessel.
# #
# #  Struttura pronta per k1 ≠ k2: basta passare griglie diverse.
# # ─────────────────────────────────────────────────────────────────────────────

function precompute_bessel_matrices(ℓ::Number,
                                    chi::AbstractVector,
                                    k1::AbstractVector,
                                    k2::AbstractVector)
    Nz = length(chi)
    N1 = length(k1)
    #N2 = length(k2)

    Bessel1 = Matrix{Float64}(undef, Nz, N1)
    #Bessel2 = Matrix{Float64}(undef, Nz, N2)

    Threads.@threads for i in 1:Nz
        @views Bessel1[i, :] = SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k1)
        #@views Bessel2[i, :] = SpecialFunctions.sphericalbesselj.(ℓ, chi[i] * k2)
    end

    #return Bessel1, Bessel2
    return Bessel1
end

# # ─────────────────────────────────────────────────────────────────────────────
# #  bessel_cheb_eval_beyond — ora delega a build_chebyshev_matrix e
# #  precompute_bessel_matrices. Argomento z_range rimosso dalla firma
# #  interna (non era usato), chi calcolato su z_cheb come nell'originale.
# # ─────────────────────────────────────────────────────────────────────────────

function bessel_cheb_eval_beyond_opt(ℓ::Number, zmin::Number, zmax::Number,
                                  kmin::Number, kmax::Number,
                                  z_range::AbstractArray, n_cheb::Int,
                                  N::Number, chi_of_z::Any)
    k  = get_clencurt_grid_z(kmin, kmax, N)
    z  = get_clencurt_grid_z(zmin, zmax, N)

    z_cheb = chebpoints(n_cheb, zmin, zmax)
    chi    = chi_of_z.(z_cheb)

    # OPT 3+4: T ora (N, n_cheb+1), costruito senza deepcopy
    T = build_chebyshev_matrix(z, zmin, zmax, n_cheb)

    # Bessel su z_cheb × k (griglia singola, come nell'originale)
    Bessel, _ = precompute_bessel_matrices(ℓ, chi, k, k)

    return T, Bessel
end

# # ─────────────────────────────────────────────────────────────────────────────
# #  compute_W_tilde — versione ottimizzata
# #
# #  OPT 1: Bessel2 precomputata fuori dal loop su p (era il collo principale).
# #  OPT 2: il loop interno è riscritto come prodotto matrice-vettore BLAS.
# #         Per ogni (p, i) fisso: T_tilde[i, p, :] = Tᵀ * (Bessel1[i,:] .* Bessel2[p,:] .* α)
# #         In forma matriciale su tutti gli i: T_tilde[:, p, :] = B_p * T
# #         dove B_p[i, k] = Bessel1[i,k] * Bessel2[p,k] * α[k].
# #  OPT 4: T è (N, n_cheb+1), Tᵀ è (n_cheb+1, N) — mul! usa dgemm BLAS.
# #  Output: T_tilde (Nz, Nz, n_cheb+1) — dimensione 1 iniziale rimossa.
# # ─────────────────────────────────────────────────────────────────────────────

function compute_W_tilde_opt(ℓ::Number, zmin::Number, zmax::Number,
                          kmin::Number, kmax::Number,
                          z_range::AbstractArray, n_cheb::Int,
                          N::Number, chi_of_z::Any)

    if zmin >= zmax
        throw(DomainError("The integration range is unphysical. Make sure zmin < zmax."))
    end

    Nz  = length(z_range)
    chi = chi_of_z.(z_range)

    # k1 e k2 sono la stessa griglia ora, pronte per essere separate in futuro
    k1 = get_clencurt_grid_z(kmin, kmax, N)
    k2 = k1   # ← sostituire con get_clencurt_grid_z(kmin2, kmax2, N2) quando serve
    α  = get_clencurt_weights_z(kmin, kmax, N)

    # T: (N, n_cheb+1) — layout colonna-major, cache-friendly per mul!
    T, _ = bessel_cheb_eval_beyond_opt(ℓ, zmin, zmax, kmin, kmax, z_range, n_cheb, N, chi_of_z)

    # OPT 1: entrambe le matrici calcolate una volta sola, O(Nz × N)
    #Bessel1, Bessel2 = precompute_bessel_matrices(ℓ, chi, k1, k2)
    Bessel1 = precompute_bessel_matrices(ℓ, chi, k1, k2)
    T_tilde = zeros(Nz, Nz, n_cheb + 1)

    # Buffer temporaneo per BLAS: (Nz, N), riusato ad ogni iterazione su p
    B_p   = Matrix{Float64}(undef, Nz, N)
    # Risultato del prodotto per uno slice p: (Nz, n_cheb+1)
    res_p = Matrix{Float64}(undef, Nz, n_cheb + 1)

    # Tᵀ: (n_cheb+1, N) — pre-trasposta per mul! (evita trasposizione interna BLAS)
    Tt = collect(T')   # (n_cheb+1, N)

    for p in 1:Nz
        # B_p[i, k] = Bessel1[i,k] * Bessel2[p,k] * α[k]
        # Bessel2[p,:] e α sono vettori di lunghezza N — broadcast su righe di B_p
        #@. B_p = Bessel1 * Bessel2[p:p, :] * α'
        @. B_p = Bessel1 * Bessel1[p:p, :] * α'

        # OPT 2: T_tilde[:, p, :] = B_p * Tᵀᵀ = B_p * T  →  (Nz,N)*(N,n_cheb+1)
        # Con Tᵀ pre-trasposto: res_p = B_p * T, shape (Nz, n_cheb+1)
        mul!(res_p, B_p, T)  # dgemm BLAS

        T_tilde[:, p, :] .= res_p
    end

    return T_tilde
end