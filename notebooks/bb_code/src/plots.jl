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
