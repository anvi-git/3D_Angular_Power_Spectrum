using Base.Threads
using LinearAlgebra

function W_tilde_computation(ℓ::Number, xmin::Number, xmax::Number, kmin::Number, kmax::Number,
                             Nk::Int, Nkp::Int, n_cheb::Int, N::Number, k_grid::AbstractVector, kp_grid::AbstractVector, x::AbstractVector)

    if xmin >= xmax
        throw(DomainError("The integration range is unphysical. Make sure xmin < xmax."))
    end

    kp = kp_grid
    #kp = get_clencurt_grid(kmin, kmax, Nkp)
    w  = get_clencurt_weights(xmin, xmax, N)

    T, Bessel1 = bessel_computation(ℓ, xmin, xmax, Nk, n_cheb, N, k_grid)
    T_tilde = zeros(eltype(w), Nk, Nkp, n_cheb + 1, 1)
    _, Bessel2 = bessel_computation(ℓ, xmin, xmax, Nkp, n_cheb, N, kp_grid)

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

function bessel_computation(ℓ::Number, xmin::Number, xmax::Number, 
                                 Nk::Int, n_cheb::Int, N::Integer, k_grid::AbstractVector)
    
    k = k_grid
    x_grid = get_clencurt_grid(xmin, xmax, N)

    T = zeros(n_cheb + 1, N)
    xx = @. (2 * x_grid - (xmax + xmin)) / (xmax - xmin)
    T[1, :] .= 1.0
    if n_cheb >= 1
        T[2, :] .= xx
    end
    for i in 3:(n_cheb + 1)
        @. T[i, :] = 2 * xx * T[i-1, :] - T[i-2, :]
    end
    
    Bessel = zeros(Nk, N)

    Threads.@threads for j in 1:Nk
        kj = k[j]
        for i in 1:N
            Bessel[j, i] = @views SpecialFunctions.sphericalbesselj.(ℓ, x_grid[i] * kj) 
        end
    end

    return T, Bessel

end