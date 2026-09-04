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

function make_k_grids(kmin, kmax, Nk, Nkp, Nkpp; sorting=false, output_dir=nothing)   

    # Log di inizio operazione (se output_dir è fornito)
    if !isnothing(output_dir)
        append_to_log(output_dir, "Generating k_grid, kp_grid, kpp_grid.")
        append_to_log(output_dir, "Parameters: kmin=$kmin, kmax=$kmax, Nk=$Nk, Nkp=$Nkp, Nkpp=$Nkpp, sorting=$sorting")
    end

    k_grid   = 10 .^ range(log10(kmin), log10(kmax), length=Nk)
    kp_grid  = 10 .^ range(log10(kmin), log10(kmax), length=Nkp)
    kpp_grid = 10 .^ range(log10(kmin), log10(kmax), length=Nkpp)

    if !isnothing(output_dir)
        append_to_log(output_dir, "grids created successfully.")
        append_to_log(output_dir, " -> k_grid:   min=$(minimum(k_grid)), max=$(maximum(k_grid)), len=$(length(k_grid))")
        append_to_log(output_dir, " -> kp_grid:  min=$(minimum(kp_grid)), max=$(maximum(kp_grid)), len=$(length(kp_grid))")
        append_to_log(output_dir, " -> kpp_grid: min=$(minimum(kpp_grid)), max=$(maximum(kpp_grid)), len=$(length(kpp_grid))")
    end
    
    return (
        k_grid   = k_grid,
        kp_grid  = kp_grid,
        kpp_grid = kpp_grid,
        sorting  = sorting
    )
end

function compute_W(ℓ::Number, xmin::Number, xmax::Number,
                             Nk::Int, Nkp::Int, n_cheb::Int, N::Int, k_grid::AbstractVector, kp_grid::AbstractVector, sorting::Bool)
    if xmin >= xmax
        throw(DomainError("The integration range is unphysical. Make sure xmin < xmax."))
    end

    if sorting == true
        w  = reverse(get_clencurt_weights(xmin, xmax, N))
    else 
        w  = get_clencurt_weights(xmin, xmax, N)
    end

    T, Bessel1 = compute_jl(ℓ, xmin, xmax, Nk, n_cheb, N, k_grid, sorting)
    T_tilde = zeros(eltype(w), Nk, Nkp, n_cheb, 1)
    #_, Bessel2 = compute_jl(ℓ, xmin, xmax, Nkp, n_cheb, N, kp_grid, sorting)

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
    #@views T_tilde[:, :, ic, 1] = Dim_Integrated .* kp_grid'
    end

    return T_tilde
end

function compute_jl(ℓ::Number, xmin::Number, xmax::Number, 
                                 Nk::Int, n_cheb::Int, N::Integer, k_grid::AbstractVector, sorting::Bool)
                             
    if sorting == true
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

function analyze_W_cheb(grid_data, x_cheb, W_x, W_cheb, c_cheb, grids, bb, paths, plot_theme)

    function cheb_eval_partial(coeffs, x, xmin, xmax, ntrunc)
        t = @. (2*x - (xmax + xmin)) / (xmax - xmin)
        T0 = ones(length(x))
        if ntrunc == 0
            return coeffs[1] .* T0
        end
        T1 = t
        S = coeffs[1] .* T0 .+ coeffs[2] .* T1
        for n in 2:ntrunc
            T2 = @. 2*t*T1 - T0
            S .+= coeffs[n+1] .* T2
            T0, T1 = T1, T2
        end
        return S
    end

    w_dir = joinpath(paths.plot_subdir, "W")
    isdir(w_dir) || mkpath(w_dir)

    # Plot 1: W(x) vs W_cheb
    p_wx = plot(grid_data.x, W_x, label = L"W(\chi)", lw = 1, ls = :dot; markercolor = :green, plot_theme.shared_style...)
    plot!(p_wx, x_cheb, W_cheb, label = L"W(\chi) \; \mathrm{on} \; \mathrm{Chebyshev} \; \mathrm{nodes}", titlefontsize = 15,
          lw = 1, ls = :dot; markercolor = :red, plot_theme.shared_style...)
    savefig(p_wx, joinpath(w_dir, "W_x_and_W_cheb.png"))

    # Calcolo W su nodi di Chebyshev (rimossa la chiamata duplicata)
    z_cheb = grid_data.z_of_χ.(x_cheb)
    W_x_on_cheb, _, _, _, _ = bb.compute_Wx(x_cheb, z_cheb, grid_data.cosmo; output_dir = paths.output_dir)

    rel_err = (W_cheb ./ W_x_on_cheb .- 1)
    rel_err_abs = abs.(rel_err)
    rel_err_pct = 100 .* rel_err
    rel_err_pct_abs = 100 .* rel_err_abs
    max_rel_err_pct = maximum(rel_err_pct_abs)

    # Plot 2: relative percentage error of Chebyshev approximation
    p_relerr = plot(1:grid_data.n_cheb, rel_err_pct,
        label = L"\frac{W_{\mathrm{Cheb}} - W_{\mathrm{true}}}{W_{\mathrm{true}}}\,[\%]",
        lw = 1.5, marker = :circle, markersize = 2, legendposition = :topright,
        xlabel = L"\chi\,[\mathrm{Mpc}/h]", ylabel = L"\mathrm{errore\ relativo}\,[\%]",
        title = "Relative error of Chebyshev approximation", size = plot_theme.size_Cl; plot_theme.shared_style...)
    hline!(p_relerr, [0.0], ls = :dash, lw = 1, color = :black, label = false)

    # Subplot combinato p1 + p2
    p1 = plot(sort(x_cheb), W_x_on_cheb,
        label = L"\sum_{n=0}^{N_{\mathrm{cheb}}-1} c_n T_n(\chi)", ls = :dash, color = :brown2)
    plot!(p1, x_cheb, W_cheb, label = L"W(\chi)", lw = 1, ls = :dot, color = :slategrey)
    plot!(p1, dpi = 1500, legendposition = :topleft, xlabel = L"\chi [Mpc/h]", ylabel = L"W(\chi) [\mathrm{Mpc}/h]",
            size = plot_theme.size_Cl; plot_theme.shared_style...)

    p2 = plot(1:grid_data.n_cheb, rel_err_pct,
            label = L"\frac{W_{\mathrm{Cheb}} - W_{\mathrm{true}}}{W_{\mathrm{true}}} [\%]",
            marker = :circle, markerstrokecolor = :black, markercolor = :black, markersize = 2,
            ls = :solid, lw = 0.5, color = :black,
            title = "relative error [%]", titlefontsize = 10, xlabel = " Chebyshev steps")
    hline!(p2, [0.0], ls = :solid, lw = 1, color = :red, label = false)

    p_combined = plot(p1, p2, layout = @layout([a; b]), dpi = 1500, legendfontsize = 15, legendposition = :outertopright,
            xlabel = [L"\chi\,[\mathrm{Mpc}/h]" "steps"], size = plot_theme.size_Cl; plot_theme.shared_style...)
    savefig(p_combined, joinpath(w_dir, "chebcoeff_combined_galaxy.png"))

    # Plot 3: c_cheb decay
    n_idx = 0:(grid_data.n_cheb - 1)
    p_decay = plot(n_idx, abs.(c_cheb),
        yscale = :log10,
        xlabel = L"n", ylabel = L"|c_n|",
        label = L"|c_n| \ \mathrm{di} \ W(\chi)",
        title = "Chebyshev coefficients decay", titlefontsize = 15, 
        marker = :circle, markercolor = :slategrey, ls = :solid, lw = 0.5; plot_theme.shared_style...)
    p_decay = plot!(p_decay, dpi = 1500)
    savefig(p_decay, joinpath(w_dir, "chebcoeff_decay_galaxy.png"))

    # Plot 4: truncation error analysis
    truncs = 1:1:length(c_cheb)-1
    errs = Float64[]
    for ntrunc in truncs
        Wn = cheb_eval_partial(c_cheb, x_cheb, grid_data.xmin, grid_data.xmax, ntrunc)
        err = maximum(abs.(Wn .- W_cheb)) / maximum(abs.(W_cheb))
        push!(errs, err)
    end

    yticks_vals = 10.0 .^ (floor(log10(minimum(errs))):ceil(log10(maximum(errs))))

    p_trunc = plot(truncs, errs,
        yscale = :log10, marker = :circle, markercolor = :slategrey, ls = :solid, lw = 0.5,
        xlabel = L"N_{trunc}",
        ylabel = L"\mathrm{err}(N_{trunc}) = \frac{\max_i |W_{N_{trunc}} - W(\chi_i)|}{\max_i |W(\chi_i)|}",
        title  = L"\mathrm{err}(N_{trunc}) = \frac{\max_i |W_{N_{trunc}} - W(\chi_i)|}{\max_i |W(\chi_i)|}",
        legend = false, titlefontsize = 15, framestyle = :box, yticks = yticks_vals, size = plot_theme.size_Cl; plot_theme.shared_style...)
    p_trunc = plot!(p_trunc, dpi = 1500)
    savefig(p_trunc, joinpath(w_dir, "chebcoeff_truncation_error_galaxy.png"))

    return (W_x_on_cheb = W_x_on_cheb, rel_err_pct = rel_err_pct, errs = errs,
            p_wx = p_wx, p_relerr = p_relerr, p_combined = p_combined, p_decay = p_decay, p_trunc = p_trunc, max_rel_err_pct = max_rel_err_pct)
end