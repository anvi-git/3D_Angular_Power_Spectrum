function bessel_cheb_eval(ℓ::Number, chi::AbstractArray, zed::AbstractArray, kmin::Number, kmax::Number, zmin::Number, zmax::Number; n_cheb::Int, N::Number)

    n_chi = length(chi)
    k = get_clencurt_grid(kmin, kmax, N)
    z = get_clencurt_grid(zmin, zmax, N)
#    zed = DataInterpolations.AkimaInterpolation(chi, z);
    z_cheb = chebpoints(n_cheb, log10(zmin), log10(zmax)) 
    c = FastChebInterp.ChebPoly(z_cheb, SA[log10(zmin)], SA[log10(zmax)])

    W = zeros(n_cheb+1,N) 
    Threads.@threads for i in 1:n_cheb+1
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.
        W[i,:] = copy_c.(log10.(z))
    end

    Bessel = zeros(n_chi, N)
    Threads.@threads for i in 1:n_chi
            Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, zed[i] * k)
    end

    return W, Bessel

end

function compute_W̃(ℓ::Number, chi::AbstractArray, zed::AbstractArray, kmin::Number, kmax::Number, zmin::Number, zmax::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
    if zmin >= zmax 
        throw(DomainError("The integration range is unphysical. Make sure zmin < zmax.")) 
    end
    
    nz = length(zed)
    n_chi = length(chi)
    z = get_clencurt_grid(zmin, zmax, N)
    w = get_clencurt_weights(kmin, kmax, N)
    k = get_clencurt_grid(kmin, kmax, N)
    kp = k
    W, Bessel1 = bessel_cheb_eval(ℓ, chi, zed, kmin, kmax, zmin, zmax; n_cheb, N)

    W_tilde = zeros(1, nz, nz, n_cheb+1)
    
    for (p, q) in enumerate(chi)
        Bessel2 = zeros(n_chi, N)
        
        Threads.@threads for i in 1:n_chi
            Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, zed[i] * kp)
        end

        α = w
         
        @tturbo for l in 1:n_cheb+1, i in 1:nz
            Cij = zero(eltype(w))
            for k in 1:N
                Cij +=  W[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
            end
            W_tilde[1,i,p,l] = Cij
        end
    end

    return W_tilde

end

function compute_W̃_chebcoeffs(
    chi_interp,       # interpolatore χ(z)
    n_z,              # funzione n(z)
    b_z,              # funzione b(z) (bias)
    D_z,              # funzione D(z) (growth factor)
    #H_z,              # funzione H(z)
    zmin::Number,
    zmax::Number;
    n_cheb::Int,
    N::Int = 2^15 + 1
)
    # Stessi nodi e pesi usati in compute_W_tilde_modes
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    t_cc   = FastTransforms.clenshawcurtisnodes(Float64, N)
    w_cc   = FastTransforms.clenshawcurtisweights(CC_obj)

    z_cc = @. (zmax + zmin) / 2 + (zmax - zmin) / 2 * t_cc
    χ_cc = chi_interp.(z_cc)
    jac  = (zmax - zmin) / 2

    # Valuta W_tilde su tutti i nodi CC
    #W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc) * H_z(z_cc)
    W_vals = @. χ_cc^2 * n_z(z_cc) * b_z(z_cc) * D_z(z_cc)

    # Costruisci T_mat: stessa logica di compute_W_tilde_modes
    c_basis = FastChebInterp.ChebPoly(zeros(n_cheb+1), SA[-1.0], SA[1.0])
    T_mat = zeros(n_cheb+1, N)
    Threads.@threads for n in 1:n_cheb+1
        basis = deepcopy(c_basis)
        basis.coefs[n] = 1.0
        T_mat[n, :] = basis.(t_cc)
    end

    # Proiezione: c_n = jac * sum_t w_cc[t] * W_vals[t] * T_n(t_cc[t])
    c_n = zeros(n_cheb+1)
    @tturbo for n in 1:n_cheb+1
        s = zero(Float64)
        for t in 1:N
            s += w_cc[t] * W_vals[t] * T_mat[n, t]
        end
        c_n[n] = jac * s
    end

    return c_n  # vettore [n_cheb]
end

function compute_Sℓ(
    I_A::AbstractArray,
    I_B::AbstractArray,
    k_array::AbstractVector;
    Pk::AbstractVector = ones(length(k_array)),
    normalization::Real = 2/pi
)
    nk = length(k_array)

    w_k = Blast.simpson_weight_array(nk)
    Δk = nk > 1 ? (last(k_array) - first(k_array)) / (nk - 1) : one(eltype(k_array))

    integrand_k = normalization .* k_array.^2 .* Pk .* w_k .* Δk

    @tullio S[li, ik1, ik2] := integrand_k[m] *
        I_A[li, ik1, m] *
        I_B[li, ik2, m]

    return S
end