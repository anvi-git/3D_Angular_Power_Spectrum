using Base.Threads
using LinearAlgebra

## BELOW: from BLAST
"""
    get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
Return the integration points in k. They are a set of 'N' Chebyshev points rescaled between 'kmin' and 'kmax'.
"""
function get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    x = (kmax - kmin) / 2 * x .+ (kmin + kmax) / 2 
    return x
end

"""
    get_clencurt_weights(kmin::Number, kmax::Number, N::Number)
Return the set of 'N' weights needed to perform the integration with the Clenshaw-Curtis quadrature rule.
The weights are rescaled between 'kmin' and 'kmax'.  
"""
function get_clencurt_weights(Xmin::Number, Xmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (Xmax - Xmin) / 2 * w

    return w
end
## ABOVE: from BLAST

function get_clencurt_grid_log(kmin::Number, kmax::Number, N::Integer)
    u_grid = Blast.get_clencurt_grid(log(kmin), log(kmax), N)  # CC nodes mapped onto [ln kmin, ln kmax]
    return exp.(u_grid)                                        # k = e^u
end

function get_clencurt_weights_log(kmin::Number, kmax::Number, N::Integer)
    k_grid = get_clencurt_grid_log(kmin, kmax, N)
    w_u = Blast.get_clencurt_weights(log(kmin), log(kmax), N)  # CC weights for ∫ du
    return w_u .* k_grid                                        # Jacobian dk = k du
end

"""
    generate_k_grids(kmin, kmax, Nk, Nkp, Nkpp; sorting) -> NamedTuple

Genera le griglie di Clenshaw-Curtis per k, kp e kpp utilizzando i limiti e i nodi forniti.
Se `sorting=true`, le griglie vengono ordinate in modo crescente.
"""
function generate_k_grids(kmin, kmax, Nk, Nkp, Nkpp; sorting=false)    
    if sorting
        k_grid   = reverse(Blast.get_clencurt_grid(kmin, kmax, Nk))
        kp_grid  = reverse(Blast.get_clencurt_grid(kmin, kmax, Nkp))
        kpp_grid = reverse(Blast.get_clencurt_grid(kmin, kmax, Nkpp))
    else
        k_grid   = Blast.get_clencurt_grid(kmin, kmax, Nk)
        kp_grid  = Blast.get_clencurt_grid(kmin, kmax, Nkp)
        kpp_grid = Blast.get_clencurt_grid(kmin, kmax, Nkpp)   
    end
    
    return (
        k_grid   = k_grid,
        kp_grid  = kp_grid,
        kpp_grid = kpp_grid,
        sorting  = sorting
    )
end

function W_tilde_computation(ℓ::Number, xmin::Number, xmax::Number, kmin::Number, kmax::Number,
                             Nk::Int, Nkp::Int, n_cheb::Int, N::Int, k_grid::AbstractVector, kp_grid::AbstractVector, x::AbstractVector, sorting::Bool)

    if xmin >= xmax
        throw(DomainError("The integration range is unphysical. Make sure xmin < xmax."))
    end

    if sorting
        w  = reverse(get_clencurt_weights(xmin, xmax, N))
    else 
        w  = get_clencurt_weights(xmin, xmax, N)
    end

    T, Bessel1 = bessel_computation(ℓ, xmin, xmax, Nk, n_cheb, N, k_grid, sorting)
    T_tilde = zeros(eltype(w), Nk, Nkp, n_cheb, 1)
    #_, Bessel2 = bessel_computation(ℓ, xmin, xmax, Nkp, n_cheb, N, kp_grid, sorting)

    # for ic in 1:n_cheb 
    #     @tturbo for ik in 1:Nk, ikp in 1:Nkp
    #         Cij = zero(eltype(w))
    #         for iN in 1:N
    #             Cij += T[ic, iN] * Bessel1[ik, iN] * Bessel2[ikp, iN] * w[iN] * kp_grid[ikp]
    #         end
    #     T_tilde[ik, ikp, ic, 1] = Cij
    #     end
    # end

    # EVEN faster? check -> is approx up to 1e-10
    for ic in 1:n_cheb
    @views Tw = T[ic, :] .* w
    #Dim_Integrated = Bessel1 * Diagonal(Tw) * Bessel2'
    Dim_Integrated = Bessel1 * Diagonal(Tw) * Bessel1'
    @views @. T_tilde[:, :, ic, 1] = Dim_Integrated * kp_grid'
    end

    return T_tilde
end

function bessel_computation(ℓ::Number, xmin::Number, xmax::Number, 
                                 Nk::Int, n_cheb::Int, N::Integer, k_grid::AbstractVector, sorting::Bool)
                             
    if sorting
        x_grid = reverse(get_clencurt_grid(xmin, xmax, N))
    else 
        x_grid = get_clencurt_grid(xmin, xmax, N)
    end
    T = zeros(n_cheb, N)
    xx = @. (2 * x_grid - (xmax + xmin)) / (xmax - xmin)
    T[1, :] .= 1.0
    if n_cheb >= 2 #1
        T[2, :] .= xx
    end
    for i in 3:(n_cheb)
        @. T[i, :] = 2 * xx * T[i-1, :] - T[i-2, :]
    end

    Bessel = zeros(Nk, N)

    Threads.@threads for j in 1:Nk
        kj = k_grid[j]
        for i in 1:N
           Bessel[j, i] = @views SpecialFunctions.sphericalbesselj.(ℓ, x_grid[i] * kj) 
        end
    end

    return T, Bessel

end