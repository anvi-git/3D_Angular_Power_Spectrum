# plot_k_grids(k_grid::AbstractVector, kp_grid::AbstractVector, kpp_grid::AbstractVector;
#               output_dir::AbstractString, plot_style::AbstractDict = Dict())
# This function generates plots of the k grid points, including histograms and heatmaps. 
# The plots are saved in the specified directory.
# ### Parameters
# - `k_grid`: An array containing the k grid points.
# - `kp_grid`: An array containing the k' grid points.
# - `kpp_grid`: An array containing the k'' grid points.
# - `output_dir`: The directory where the plots will be saved.
# - `plot_style`: A dictionary with plot style options (e.g., line color, title, legend position).

# ### Returns
# - `hist_k`: Histogram of the k grid points.
# - `hist_kp`: Histogram of the k' grid points.
# - `hist_kpp`: Histogram of the k'' grid points.

function plot_k_grids(k_grid::AbstractVector, kp_grid::AbstractVector, kpp_grid::AbstractVector;
                      output_dir::AbstractString, plot_style::AbstractDict = Dict())

    Nk   = length(k_grid)
    Nkp  = length(kp_grid)
    Nkpp = length(kpp_grid)

    bins_k   = round(Int, sqrt(Nk))
    bins_kp  = round(Int, sqrt(Nkp))
    bins_kpp = round(Int, sqrt(Nkpp))

    plots_dir = joinpath(output_dir, "plots")
    k_grid_dir = joinpath(plots_dir, "k_grid")
    isdir(k_grid_dir) || mkpath(k_grid_dir)
    color = get(plot_style, :color, :slategrey)
    hist_k_log = histogram(log10.(k_grid), bins = bins_k,
              xlabel = L"\mathrm{log_{10}}(k \; (h/\mathrm{Mpc}))", ylabel = "Number of points", color = color,
              title = "Distribution of k grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "k_grid_distribution_log.png"))

    hist_k = histogram(k_grid, bins = bins_k,
              xlabel = L"k \; (h/\mathrm{Mpc})", ylabel = "Number of points", color = color,
              title = "Distribution of k grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "k_grid_distribution.png"))

    p_combined = plot(hist_k_log, hist_k, layout = (2, 1), size = (800, 900), dpi = 150)
    savefig(p_combined, joinpath(k_grid_dir, "k_grid_distribution_combined.png"))
    
    hist_kp = histogram(kp_grid, bins = bins_kp,
              xlabel = L"k' \; (h/\mathrm{Mpc})", ylabel = "Number of points",
              title = "Distribution of k' grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "kp_grid_distribution.png"))

    hist_kpp = histogram(kpp_grid, bins = bins_kpp,
              xlabel = L"k'' \; (h/\mathrm{Mpc})", ylabel = "Number of points",
              title = "Distribution of k'' grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "kpp_grid_distribution.png"))

    return hist_k_log, hist_k, hist_kp, hist_kpp
end

function plot_heatmaps(
    W_final_gal,
    k_grid,
    kp_grid,
    plot_theme,
    paths;
    il::Int = 1
)
    xticks_vals = 10.0 .^ (floor(Int, log10(minimum(k_grid))) :ceil(Int,  log10(maximum(k_grid))))
    yticks_vals = 10.0 .^ (floor(Int, log10(minimum(kp_grid))) :ceil(Int,  log10(maximum(kp_grid))))
    Wnorm = W_final_gal[il, :, :] / maximum(W_final_gal[il, :, :])

    p1 = heatmap(
        k_grid, kp_grid,
        Wnorm,
        title = L"W_{final}^{gg}" * "(at fixed " * L"ℓ = 2)",
        xscale = :log10,
        yscale = :log10,
        xlabel = L"log_{10}(k) \; [\mathrm{h/Mpc}]",
        ylabel = L"log_{10}(k_1) \; [\mathrm{h/Mpc}]",
        size = plot_theme.size_heatmap,
        c = plot_theme.c,
        xticks = xticks_vals,
        yticks = yticks_vals;
        plot_theme.shared_style...
    )

    savefig(p1, joinpath(paths.plot_subdir, "Sl_plots/heatmaps/heatmap_k_k1_log.png"))

    p2 = heatmap(
        k_grid, kp_grid,
        Wnorm,
        title = L"W_{final}^{gg}" * "(at fixed " * L"ℓ = 2)",
        xlabel = L"k",
        ylabel = L"k_1",
        size = plot_theme.size_heatmap,
        c = plot_theme.c;
        plot_theme.shared_style...
    )

    savefig(p2, joinpath(paths.plot_subdir, "Sl_plots/heatmaps/heatmap_k_k1.png"))

    return (; p1, p2, xticks_vals, yticks_vals)
end

function plot_theory_Pk(
    grid_data,
    plot_theme;
    pk_file = "bb_code/data/pk.npz",
    χ1::Float64 = 1000.0,
    χ2::Float64 = 1000.0
)

    pk_dict = npzread(pk_file)
    Pklin = pk_dict["pk_lin"]
    Pknonlin = pk_dict["pk_nl"]
    k_1k = pk_dict["k"]
    z_pk = pk_dict["z"]

    y_pk = LinRange(log10(first(k_1k)), log10(last(k_1k)), length(k_1k))
    x_pk = LinRange(first(z_pk), last(z_pk), length(z_pk))

    InterpPmm = Interpolations.interpolate(
        log10.(Pklin),
        BSpline(Cubic(Line(OnGrid())))
    )
    InterpPmm = scale(InterpPmm, (x_pk, y_pk))
    InterpPmm = Interpolations.extrapolate(InterpPmm, Line())

    InterpPmm_nl = Interpolations.interpolate(
        log10.(Pknonlin),
        BSpline(Cubic(Line(OnGrid())))
    )
    InterpPmm_nl = scale(InterpPmm_nl, x_pk, y_pk)
    InterpPmm_nl = Interpolations.extrapolate(InterpPmm_nl, Line())

    power_spectrum(k, χ1, χ2) = @. sqrt(
        10^InterpPmm(grid_data.z_of_χ(χ1), log10(k)) *
        10^InterpPmm(grid_data.z_of_χ(χ2), log10(k))
    )

    power_spectrum_nl(k, χ1, χ2) = @. sqrt(
        10^InterpPmm_nl(grid_data.z_of_χ(χ1), log10(k)) *
        10^InterpPmm_nl(grid_data.z_of_χ(χ2), log10(k))
    )

    pk_lin_vals = power_spectrum.(k_1k, χ1, χ2)
    pk_nl_vals  = power_spectrum_nl.(k_1k, χ1, χ2)

    p = plot(
        k_1k, pk_lin_vals,
        label = "Linear P(k)",
        xscale = :log10,
        yscale = :log10,
        title = L"$P(k)$ at fixed $\chi_1,\chi_2$",
        titlefontsize = 20,
        xlabel = L"$k \; (h/\mathrm{Mpc})$",
        ylabel = L"$P(k) \; ((\mathrm{Mpc}/h)^3)$",
        labelfontsize = 15,
        dpi = 1500
    )

    plot!(
        p,
        k_1k, pk_nl_vals,
        label = "Non-linear P(k)",
        xscale = :log10,
        yscale = :log10,
        size = plot_theme.size_Cl;
        plot_theme.shared_style...
    )

    return (; p, k_1k, pk_lin_vals, pk_nl_vals, InterpPmm, InterpPmm_nl, power_spectrum, power_spectrum_nl)
end

function plot_Sl_kp_kpp(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                         i::Int=150, j::Int=150, save_fig::Bool=true, showfig::Bool=false)
    
    k_1 = grids.k_grid[i]
    k_2 = grids.kp_grid[j]

    p = plot(grid_data.ℓ, S_lkk_gg[i, j, :],
        color = :black,
        label = L"k_1 = k_{2} = %$(round(k_1, digits=5)) \; \mathrm{h/Mpc}"
    )

    plot!(p,
        xlabel = L"\ell",
        ylabel = L"S_\ell (\mathrm{Mpc}/h)^2",
        # xscale = :log10,
        minorticks = true
    )

    plot!(p;
        label = "Beyond BLAST",
        size = plot_theme.size_Cl,
        title = L"S_{\ell}^{gg} (k_1 = k_{2} = %$(round(k_1, digits=5)) \; \mathrm{h/Mpc})",
        titlefontsize = 20,
        titleposition = :left,
        dpi = 1500,
        plot_theme.shared_style...
    )

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        savefig(p, joinpath(paths.plot_subdir, "Sl_plots/ell_vs_Sl_kp_kpp_$(i)_$(j).png"))
    end

    if showfig
        display(p)
    end

    return p
end

function plot_Sl_fixed_l_kpp(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                            il::Int=1, ikpp::Int=5, save_fig::Bool=true, showfig::Bool=false)
    
    k_2 = grids.kp_grid[ikpp]
    ell_val = grid_data.ℓ[il]

    p = plot(sort(grids.k_grid), S_lkk_gg[:, ikpp, il], 
        label = L"k_{2} = %$(round(k_2, digits=5)) \; \mathrm{h/Mpc}, \; \ell = %$(ell_val)",
        linestyle = :solid, 
        lw = 1.5
    )

    plot!(p,
        xlabel = L"k_{1} \; (h/\mathrm{Mpc})",
        ylabel = L"S_\ell \; (\mathrm{Mpc}/h)^2",
        minorticks = true,
        xscale = :log10
    )

    plot!(p,
        title = L"S_\ell^{gg} = \int dk k^2 P(k) \int \tilde W(k, k_1) \int \tilde W(k, k_2) (\mathrm{varying} \; k_1 (k_1), \mathrm{at} \; \mathrm{fixed} \; k_2 (k_{2}))",
        titlefontsize = 20,
        label = L"\ell=%$(ell_val), k_{2}=%$(round(k_2, digits=3)) \; \mathrm{h/Mpc}",
        titleposition = :left, 
        labelfontsize = 20, 
        legendposition = :outertopright, 
        size = plot_theme.size_Cl; 
        plot_theme.shared_style...
    )

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        fname = "kpp_vs_Sl_idxell_$(il)_kpp_$(round(k_2, digits=3)).png"
        savefig(p, joinpath(paths.plot_subdir, "Sl_plots", fname))
    end

    if showfig
        display(p)
    end

    return p
end

function plot_Sl_fixed_l_varying_k(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                                   il::Int=1, step::Int=20, normalize::Bool=true, 
                                   save_fig::Bool=true, showfig::Bool=false)

    ell_val = grid_data.ℓ[il]
    p = plot()

    # Pre-generate color gradient across the total grid length
    color_palette = palette(:roma, grid_data.Nkp)

    for ip in 1:step:grid_data.Nkp
        kp = grids.k_grid[ip]
        
        # Extract signal curve for current k_2 index and l index
        signal = S_lkk_gg[:, ip, il]
        y_data = normalize ? signal ./ maximum(signal) : signal

        plot!(p, sort(grids.k_grid), y_data, 
            color = color_palette[ip],
            label = L"k_1 = k_{2} = %$(round(kp, digits=5)) \; \mathrm{h/Mpc}",
            linestyle = :solid,
            lw = 1.5
        )
    end

    plot!(p,
        xlabel = L"k_{1} \; (h/\mathrm{Mpc})",
        ylabel = normalize ? L"S_\ell / \max(S_\ell)" :  L"S_\ell (\mathrm{Mpc}/h)^2",
        #xscale = :log10,
        minorticks = true
    )

    plot!(p,
        title = L"S_\ell^{gg} \; \mathrm{for} \; \ell = %$(ell_val)",
        legendposition = :outertopright, 
        size = plot_theme.size_Cl, 
        titlefontsize = 20, 
        titleposition = :left; 
        plot_theme.shared_style...
    )
    if step < 10
        plot!(p, legend = false)
        scatter!(p, [NaN], [NaN], 
        zcolor = [minimum(grids.k_grid)], 
        clims = (minimum(grids.k_grid), maximum(grids.k_grid)),
        c = :viridis,
        colorbar_title = L"k_1 \; (\mathrm{h/Mpc})",
        label = ""
        )
    end

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        suffix = normalize ? "_normalized" : ""
        fname = "kp_vs_Sl_idxell_$(il)_kpp_fixed$(suffix).png"
        savefig(p, joinpath(paths.plot_subdir, "Sl_plots", fname))
    end

    if showfig
        display(p)
    end

    return p
end

function animate_Sl_fixed_l_varying_k(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                                       il::Int=1, step::Int=20, normalize::Bool=false, 
                                       fps::Int=1, save_fig::Bool=true, showfig::Bool=false)

    ell_val = grid_data.ℓ[il]
    
    plt = plot(
        xlabel = L"k_{1} \; (h/\mathrm{Mpc})",
        ylabel = normalize ? L"S_\ell / \max(S_\ell)" : L"S_\ell (\mathrm{Mpc}/h)^2",
        xscale = :log10,
        minorticks = true,
        title = L"S_\ell^{gg} \; \mathrm{for} \; \ell = %$(ell_val)",
        legendposition = :outertopright, 
        size = plot_theme.size_Cl, 
        titlefontsize = 20, 
        titleposition = :left; 
        plot_theme.shared_style...
    )

    color_palette = cgrad(:seaborn_icefire_gradient, grid_data.Nkp)

    anim = @animate for ip in 1:step:grid_data.Nkp
        kp = grids.k_grid[ip]
        signal = S_lkk_gg[:, ip, il]
        y_data = normalize ? signal ./ maximum(signal) : signal

        plot!(
            plt,
            grids.k_grid, 
            y_data, 
            color = color_palette[ip],
            label = L"k_1 = k_{2} = %$(round(kp, digits=5)) \; \mathrm{h/Mpc}",
            linestyle = :solid, 
            lw = 1.5
        )
    end

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        suffix = normalize ? "_normalized" : ""
        gif_path = joinpath(paths.plot_subdir, "Sl_plots", "kp_vs_Sl_idxell_$(il)_evolution$(suffix).gif")
        gif(anim, gif_path, fps = fps)
    end

    if showfig
        display(anim)
    end

    return anim
end

function plot_Sl_fixed_kpp_varying_l(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                                     ipp::Int=4, step::Int=20, normalize::Bool=false, 
                                     save_fig::Bool=true, showfig::Bool=false)

    kpp = grids.kp_grid[ipp]
    p = plot()

    n_ell = length(grid_data.ℓ)
    color_palette = cgrad(:seaborn_icefire_gradient, n_ell)

    for il in 1:step:n_ell
        signal = S_lkk_gg[:, ipp, il]
        y_data = normalize ? signal ./ maximum(signal) : signal

        plot!(p, grids.kp_grid, y_data,
            color = color_palette[il],
            label = L"k_{2} = %$(round(kpp, digits=5)) \; \mathrm{h/Mpc}, \; \ell = %$(grid_data.ℓ[il])",
            linestyle = :solid, 
            lw = 1.5
        )
    end

    plot!(p,
        xlabel = L"k_{1} \; (h/\mathrm{Mpc})",
        ylabel = normalize ? L"S_\ell / \max(S_\ell)" : L"S_\ell (\mathrm{Mpc}/h)^2",
        xscale = :log10,
        minorticks = true
    )

    plot!(p,
        title = L"S_\ell^{gg} = \int dk k^2 P(k) \int \tilde W(k, k_1) \int \tilde W(k, k_2) (\mathrm{at} \; \mathrm{different} \; \ell)",
        titlefontsize = 20,
        titleposition = :left,
        legendposition = :outertopright, 
        size = plot_theme.size_Cl; 
        plot_theme.shared_style...
    )

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        suffix = normalize ? "_normalized" : ""
        fname = "kp_vs_Sl_idxkpp_$(ipp)_varying_l$(suffix).png"
        savefig(p, joinpath(paths.plot_subdir, "Sl_plots", fname))
    end

    if showfig
        display(p)
    end

    return p
end

using Plots

function animate_Sl_fixed_kpp_varying_l(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                                        ipp::Int=140, step::Int=30, normalize::Bool=false, 
                                        fps::Int=1, save_fig::Bool=true, showfig::Bool=false)

    kpp = grids.kp_grid[ipp]
    n_ell = length(grid_data.ℓ)
    
    # Initialize the base plot outside the loop
    plt = plot(
        xlabel = L"k_{1} \; (h/\mathrm{Mpc})",
        ylabel = normalize ? L"S_\ell / \max(S_\ell)" : L"S_\ell (\mathrm{Mpc}/h)^2",
        xscale = :log10,
        minorticks = true,
        title = L"S_\ell^{gg} = \int dk k^2 P(k) \int \tilde W(k, k_1) \int \tilde W(k, k_2) \; (\mathrm{at} \; \mathrm{different} \; \ell)",
        titlefontsize = 16,
        titleposition = :left, 
        legendposition = :outertopright, 
        size = plot_theme.size_Cl; 
        plot_theme.shared_style...
    )

    color_palette = cgrad(:seaborn_icefire_gradient, n_ell) 

    # Generate the animation frame by frame
    anim = @animate for il in 1:step:n_ell
        signal = S_lkk_gg[:, ipp, il]
        y_data = normalize ? signal ./ maximum(signal) : signal

        plot!(
            plt,
            grids.kp_grid, 
            y_data,
            color = color_palette[il],
            label = L"k_{2} = %$(round(kpp, digits=5)) \; \mathrm{h/Mpc}, \; \ell = %$(grid_data.ℓ[il])",
            linestyle = :solid
        )
    end

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        suffix = normalize ? "_normalized" : ""
        # Cleaned up filename to avoid excessively long float strings
        fname = "kpp_vs_Sl_idxell_loop_on_ell_kp_$(round(kpp, digits=3))$(suffix).gif"
        gif_path = joinpath(paths.plot_subdir, "Sl_plots", fname)
        
        gif(anim, gif_path, fps = fps)
    end
    
    if showfig
        display(anim)
    end
    
    return anim
end

function plot_l_lplus1_Sl(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                          ip::Int=150, ipp::Int=150, save_fig::Bool=true, showfig::Bool=false)
    
    k_1 = grids.k_grid[ip]
    k_2 = grids.kp_grid[ipp]

    # Calculate \ell(\ell+1) S_\ell
    y_data = S_lkk_gg[ip, ipp, :] .* grid_data.ℓ .* (grid_data.ℓ .+ 1)

    p = plot(
        grid_data.ℓ,
        y_data,
        xlabel = L"\ell",
        ylabel = L"\ell(\ell+1)S_\ell (\mathrm{Mpc}/h)^2",
        label = "Beyond BLAST",
        legend = false,
        colorbar = true,
        colorbar_title = L"j",
        clims = (1, grid_data.Nkp),
        title = L"\ell(\ell+1)S_\ell^{gg}",
        size = plot_theme.size_Cl,
        titlefontsize = 20,
        xscale = :log10,
        titleposition = :left;
        plot_theme.shared_style...
    )

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        fname = "l_lplus1_Sl_ip_$(ip)_ipp_$(ipp).png"
        savefig(p, joinpath(paths.plot_subdir, "Sl_plots", fname))
    end

    if showfig
        display(p)
    end

    return p
end

function plot_Sl_diagonal_ell_evolution(grids, S_lkk_gg, grid_data, paths, plot_theme; 
                                       nshow::Int=100, normalize::Bool=true, 
                                       save_fig::Bool=true, showfig::Bool=false)

    # Extract diagonals across all ell indices: shape (Nkp, N_ell)
    n_ell = size(S_lkk_gg, 3)
    diagS = [diag(S_lkk_gg[:, :, i]) for i in 1:n_ell]

    # Select sampled indices across the ell range
    idx = round.(Int, range(1, n_ell, length=min(nshow, n_ell)))

    p = plot(
        xlabel = L"k \; (h/\mathrm{Mpc})",
        ylabel = normalize ? L"S_\ell(k,k) / \max(S_\ell(k,k))" : L"S_\ell(k,k) (\mathrm{Mpc}/h)^2",
        title  = L"S_\ell^{gg}(k,k)\ \mathrm{for\ different}\ \ell\ \mathrm{values}",
        colorbar_title = L"\ell",
        clims = (grid_data.ℓ[1], grid_data.ℓ[end]),
        legend = false,
        colorbar = true,
        size = plot_theme.size_Cl; 
        plot_theme.shared_style...
    )

    # Palette mapped to the number of sampled lines
    colors = cgrad(:seaborn_icefire_gradient, length(idx))

    for (k, ii) in enumerate(idx)
        curve = diagS[ii]
        y_data = normalize ? curve ./ maximum(curve) : curve

        plot!(p, grids.kp_grid, y_data, 
            line_z = grid_data.ℓ[ii],
            color = colors[k],
            lw = 1.5,
            label = L"\ell = %$(grid_data.ℓ[ii])"
        )
    end

    if save_fig
        mkpath(joinpath(paths.plot_subdir, "Sl_plots"))
        suffix = normalize ? "_normalized" : ""
        fname = "Sl_diagonal_ell_evolution$(suffix).png"
        savefig(p, joinpath(paths.plot_subdir, "Sl_plots", fname))
    end

    if showfig
        display(p)
    end

    return p
end