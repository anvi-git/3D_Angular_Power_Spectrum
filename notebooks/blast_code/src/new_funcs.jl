function compute_W_tilde_modes(ℓ::Number,
                               k_arr::AbstractVector,
                               k_other_arr::AbstractVector,
                               z_nodes::AbstractVector,
                               chi_interp,
                               zmin::Number,
                               zmax::Number;
                               n_cheb::Int = length(z_nodes),
                               N::Int = 2^(15) + 1)

    nk  = length(k_arr)
    nk2 = length(k_other_arr)

    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    t_cc   = FastTransforms.clenshawcurtisnodes(Float64, N)
    w_cc   = FastTransforms.clenshawcurtisweights(CC_obj)

    z_cc = @. (zmax + zmin) / 2 + (zmax - zmin) / 2 * t_cc
    χ_cc = chi_interp.(z_cc)

    c_basis = FastChebInterp.ChebPoly(zeros(n_cheb), SA[-1.0], SA[1.0])
    T_mat = zeros(n_cheb, N)

    Threads.@threads for n in 1:n_cheb
        basis = deepcopy(c_basis)
        basis.coefs .*= 0
        basis.coefs[n] = 1.0
        T_mat[n, :] = basis.(t_cc)
    end

    Bk  = zeros(nk, N)
    Bok = zeros(nk2, N)

    Threads.@threads for i in 1:nk
        Bk[i, :] = SpecialFunctions.sphericalbesselj.(ℓ, k_arr[i] .* χ_cc)
    end

    Threads.@threads for i in 1:nk2
        Bok[i, :] = SpecialFunctions.sphericalbesselj.(ℓ, k_other_arr[i] .* χ_cc)
    end

    jac = (zmax - zmin) / 2
    W_basis = zeros(Float64, nk2, nk, n_cheb)

    @tturbo for n in 1:n_cheb, ik2 in 1:nk2, ik in 1:nk
        s = zero(Float64)
        for t in 1:N
            s += w_cc[t] * T_mat[n, t] * Bok[ik2, t] * Bk[ik, t]
        end
        W_basis[ik2, ik, n] = jac * s
    end

    return W_basis
end