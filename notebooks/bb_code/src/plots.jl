function plot_k_grids(k_grid::AbstractVector, kp_grid::AbstractVector, kpp_grid::AbstractVector;
                      output_dir::AbstractString, plot_style::AbstractDict = Dict())

    Nk   = length(k_grid)
    Nkp  = length(kp_grid)
    Nkpp = length(kpp_grid)

    bins_k   = round(Int, sqrt(Nk))
    bins_kp  = round(Int, sqrt(Nkp))
    bins_kpp = round(Int, sqrt(Nkpp))

    # Assicura che la cartella di output esista
    plots_dir = joinpath(output_dir, "plots")
    k_grid_dir = joinpath(plots_dir, "k_grid")
    isdir(k_grid_dir) || mkpath(k_grid_dir)

    hist_k = histogram(k_grid, bins = bins_k,
              xlabel = L"k \; (h/\mathrm{Mpc})", ylabel = "Number of points", 
              title = "Distribution of k grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "k_grid_distribution.png"))

    hist_kp = histogram(kp_grid, bins = bins_kp,
              xlabel = L"k' \; (h/\mathrm{Mpc})", ylabel = "Number of points",
              title = "Distribution of k' grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "kp_grid_distribution.png"))

    hist_kpp = histogram(kpp_grid, bins = bins_kpp,
              xlabel = L"k'' \; (h/\mathrm{Mpc})", ylabel = "Number of points",
              title = "Distribution of k'' grid points", legend = false; plot_style...)
    savefig(joinpath(k_grid_dir, "kpp_grid_distribution.png"))

    return hist_k, hist_kp, hist_kpp
end

function plot_heatmaps(
    W_final_gal,
    k_grid,
    kp_grid,
    plot_theme,
    paths;
    il::Int = 1
)

    xticks_vals = 10.0 .^ (
        floor(Int, log10(minimum(k_grid))) :
        ceil(Int,  log10(maximum(k_grid)))
    )

    yticks_vals = 10.0 .^ (
        floor(Int, log10(minimum(kp_grid))) :
        ceil(Int,  log10(maximum(kp_grid)))
    )

    Wnorm = W_final_gal[il, :, :] / maximum(W_final_gal[il, :, :])

    p1 = heatmap(
        k_grid, kp_grid,
        Wnorm,
        title = L"W_{final}^{gg}" * "(at fixed " * L"ℓ)",
        xscale = :log10,
        yscale = :log10,
        xlabel = L"log_{10}(k) \; [\mathrm{h/Mpc}]",
        ylabel = L"log_{10}(k_p) \; [\mathrm{h/Mpc}]",
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
        title = L"W_{final}^{gg}" * "(at fixed " * L"ℓ)",
        xlabel = L"k",
        ylabel = L"k_p",
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
    k_pk = pk_dict["k"]
    z_pk = pk_dict["z"]

    y_pk = LinRange(log10(first(k_pk)), log10(last(k_pk)), length(k_pk))
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

    pk_lin_vals = power_spectrum.(k_pk, χ1, χ2)
    pk_nl_vals  = power_spectrum_nl.(k_pk, χ1, χ2)

    p = plot(
        k_pk, pk_lin_vals,
        label = "Linear P(k)",
        xscale = :log10,
        yscale = :log10,
        title = L"$P(k)$ at fixed $\chi_1,\chi_2$",
        titlefontsize = 20,
        xlabel = L"$k \; (h/\mathrm{Mpc})$",
        ylabel = L"$P(k) \; ((\mathrm{Mpc}/h)^3)$",
        labelfontsize = 15
    )

    plot!(
        p,
        k_pk, pk_nl_vals,
        label = "Non-linear P(k)",
        xscale = :log10,
        yscale = :log10,
        size = plot_theme.size_Cl;
        plot_theme.shared_style...
    )

    return (; p, k_pk, pk_lin_vals, pk_nl_vals, InterpPmm, InterpPmm_nl, power_spectrum, power_spectrum_nl)
end
