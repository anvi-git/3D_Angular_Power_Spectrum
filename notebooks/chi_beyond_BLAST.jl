ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0;
using Pkg
Pkg.activate("blast_code"; io=devnull)
Pkg.resolve(; io=devnull)
Pkg.instantiate(; io=devnull)

using Revise 

using Base.Threads, NPZ, DataInterpolations, Interpolations, FastChebInterp
using BenchmarkTools, FFTW, FastTransforms, Dates, TOML, Plots, Plots.Measures
using QuadGK, LaTeXStrings, Tullio, StaticArrays, LoopVectorization, LinearAlgebra
using Unitful, SpecialFunctions, DifferentialEquations, Cosmology, NumericalIntegration
using CSV, DataFrames, JSON, OrderedCollections
;

include("blast_code/src/Blast.jl")
include("blast_code/src/blast_tutorials.jl")
using .Blast
using .blast_tutorials
include("blast_code/src/galaxy_galaxy.jl")
include("blast_code/src/shear_shear.jl")
using .galaxy_galaxy
using .shear_shear

include("blast_code/src/Blast.jl")
include("blast_code/src/blast_tutorials.jl")
using .Blast
using .blast_tutorials
include("blast_code/src/galaxy_galaxy.jl")
include("blast_code/src/shear_shear.jl")
using .galaxy_galaxy
using .shear_shear

timestamp = Dates.format(now(), "yyyy_mm_dd_HHMMSS")
output_dir = "out/runs/run_$timestamp"
plot_subdir = joinpath(output_dir, "plots")
Sl_plots = joinpath(plot_subdir, "Sl_plots")
Kernel_plots = joinpath(plot_subdir, "kernels")
quantity_subdir = joinpath(output_dir, "quantities")
chebcoefs = joinpath(quantity_subdir, "chebcoefs")
Sl = joinpath(quantity_subdir, "Sl")
mkpath(output_dir)
mkpath(plot_subdir)
mkpath(quantity_subdir)
mkpath(chebcoefs)
mkpath(Sl)
mkpath(Sl_plots)
mkpath(Kernel_plots)
println("Folders in place: ", output_dir, ", ", plot_subdir, " and ", quantity_subdir)

#Background quantities
z_b = npzread("blast_code/data/background/z.npy") # array 
x_b = npzread("blast_code/data/background/chi.npy") # array
n5k_bins = npzread("blast_code/data/dNdzs_fullwidth.npz") # array
# using Akima interpolation
z_of_χ = DataInterpolations.AkimaInterpolation(z_b, x_b); # z(χ)
chi_of_z = DataInterpolations.AkimaInterpolation(x_b, z_b); # χ(z)

cosmo = Blast.FlatΛCDM()
N = 2^15+1                           # Number of comoving distance points
xmin = 26                            # minimum comoving distance [Mpc/h]
xmax = 7000                          # maximum comoving distance [Mpc/h]
x = LinRange(xmin, xmax, N)          # in Mpc/h: from 26 Mpc/h to 7000 Mpc/h
z = z_of_χ.(x)
zmin = minimum(z)                    # minimum redshift
zmax = maximum(z)                    # maximum redshift
ℓ = LinRange(2, 200, 100)
Nk = 2^10
Nkp = 2^10
Nkpp = 2^10
kmax = 200/13                        # maximum wavenumber (small scales) [h/Mpc]
kmin = 2.5/xmax                      # minimum wavenumber (large scales) [h/Mpc]
n_cheb = 200                         # number of Chebyshev nodes 
;

size_plot = (1200, 600)
size_heatmap = (600, 600)
dpi = 300
leftmargin = 5Plots.mm
bottommargin = 5Plots.mm
rightmargin = 5Plots.mm
topmargin = 5Plots.mm
titlefontsize = 15
yguidefontsize =15 
xguidefontsize =15
legendfontsize = 10
c = :curl
colors = palette(:batlow, Nkp)
;

# b1 = 0.98 .+ 1.24 .* z .- 1.72 .* z.^2 .+ 1.28 .* z.^3
# b2 = 1.0 .+ sqrt.(1.0 .+ z)
# plot(z, b1, label = "b1", xlabel = L"z", ylabel = L"bias", title = "Bias as a function of redshift", legendfontsize = legendfontsize, dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin)
# plot!(z, b2, label = "b2")

sorting = false
if sorting == true
    k_grid_sort = sort(Blast.get_clencurt_grid(kmin, kmax, Nk))
    kp_grid_sort = sort(Blast.get_clencurt_grid(kmin, kmax, Nkp))
    kpp_grid_sort = sort(Blast.get_clencurt_grid(kmin, kmax, Nkpp))
    k_grid = k_grid_sort
    kp_grid = kp_grid_sort
    kpp_grid = kpp_grid_sort
else
    k_grid = Blast.get_clencurt_grid(kmin, kmax, Nk)
    kp_grid = Blast.get_clencurt_grid(kmin, kmax, Nkp)
    kpp_grid = Blast.get_clencurt_grid(kmin, kmax, Nkpp)
end
;

# plot(k_grid, kp_grid, kpp_grid, 
#     label = ["k" "k'" "k''"], 
#     xscale = :log10, yscale = :log10, zscale = :log10,
#     size = size_plot,
#     xaxis = L"k [h/Mpc]", yaxis = L"k' [h/Mpc]", zaxis = L"k'' [h/Mpc]",
#     legendfontsize = legendfontsize, dpi = dpi, 
#     leftmargin = leftmargin, bottommargin = bottommargin)

using OrderedCollections

params_run = OrderedDict(
    "number of points for the integration of the integrals of W, N" => N,
    "minimum comoving distance xmin" => xmin,
    "maximum comoving distance xmax" => xmax,
    "minimum redshift zmin" => zmin,
    "maximum redshift zmax" => zmax,
    "minimum wavenumber kmin" => kmin,
    "maximum wavenumber kmax" => kmax,
    "number of Chebyshev nodes n_cheb" => n_cheb,
    "number of ℓ values" => length(ℓ),
    "number of k points Nk" => Nk,
    "number of kp points Nkp" => Nkp,
    "number of kpp points Nkpp" => Nkpp
)

mkpath(output_dir)

open(joinpath(output_dir, "config.txt"), "w") do dictio
    println(dictio, "=== Run parameters ===")
    for (key, val) in params_run
        println(dictio, key, " => ", val)
    end
    println(dictio)
    println(dictio, "=== Additional info ===")
    println(dictio, "comoving distance array x has size: ", length(x))
    println(dictio, "redshift array z has size: ", length(z))

    if k_grid != kp_grid
        println(dictio, "k_grid and kp_grid are different")
    else
        println(dictio, "k_grid and kp_grid are the same")
    end
    if k_grid != kpp_grid
        println(dictio, "k_grid and kpp_grid are different")
    else
        println(dictio, "k_grid and kpp_grid are the same")
    end
    if kp_grid != kpp_grid
        println(dictio, "kp_grid and kpp_grid are different")
    else
        println(dictio, "kp_grid and kpp_grid are the same")
    end

    if sorting == false
        println(dictio, "k_grid is not sorted")
    else
        println(dictio, "k_grid is sorted")
    end
end

println("Data brought to safety in folder: ", output_dir)


#shear_prefac_W = shear_shear.shear_prefactor(x, z, cosmo; output_dir=output_dir);
#shear_prefac_W_cheb = shear_shear.shear_prefactor_cheb(zmin, zmax, n_cheb, z, shear_prefac_W; output_dir=output_dir);
#cheb_coeff_shear = shear_shear.compute_prefactor_chebcoeffs(shear_prefac_W_cheb; output_dir=output_dir);

gal_prefact_W, bias, growth, Hubble_param, nz_norm = galaxy_galaxy.galaxy_prefactor(x, z, cosmo; output_dir=output_dir);
gal_prefact_W_cheb = galaxy_galaxy.galaxy_prefactor_cheb(xmin, xmax, n_cheb, z, x, bias, growth, Hubble_param, nz_norm; output_dir=output_dir);
cheb_coeff_gal = galaxy_galaxy.compute_prefactor_chebcoeffs(gal_prefact_W_cheb; output_dir=output_dir);

#println("Size of cheb_coeffs must be the same for shear and galaxy")
#println("size of cheb_coeff_shear: ", size(cheb_coeff_shear))
println("size of cheb_coeff_gal: ", size(cheb_coeff_gal))

# W_tilde = zeros(Nk, Nkp, n_cheb, length(ℓ))
# W_tilde = npzread("/Users/anvi/Desktop/cosmo/notebooks/out/W_tilde.npy")
# ;

W_tilde = zeros(Nk, Nkp, n_cheb, length(ℓ))
elapsed_time = zeros(length(ℓ))
println("Dimensions of W_tilde: ", size(W_tilde))
for i in eachindex(ℓ)
    # i want the elapsed time for each l 
    t_0 = time()
    W_tilde[:, :, :, i] .= Blast.W_tilde_computation(ℓ[i], zmin, zmax, kmin, kmax,
                                                 Nk, Nkp, n_cheb - 1, N, k_grid, kp_grid, x)
    t_end = time()
    elapsed_time[i] = t_end - t_0
    println("Finished computing W_tilde for ℓ = $(ℓ[i]). Time elapsed $(round(elapsed_time[i], digits=2))s")
end
;

#npzwrite(joinpath(output_dir, "W_tilde.npy"), W_tilde)

println("Size of W_tilde: ", size(W_tilde))
#Size of W_tilde: (Nk, Nkp, Ncheb, Nl)

# println("min = ", minimum(W_tilde), ", max = ", maximum(W_tilde))
# println("mean = ", mean(W_tilde), ", std = ", std(W_tilde))
# if  any(isnan, W_tilde)
#     println("There are NaN values in W_tilde")
# end
# if any(isinf, W_tilde)
#     println("There are infinite values in W_tilde")
# end
# if any(iszero, W_tilde)
#     println("There are zero values in W_tilde")
# end


#npzwrite(joinpath(quantity_subdir, "chebcoefs/cheb_coeff_gal.npy"), cheb_coeff_gal)

@tullio W_final_gal[il, ik, ikp] := W_tilde[ik, ikp, ic, il] * cheb_coeff_gal[ic];

println("Size of W_final_gal: ", size(W_final_gal))
# [ell, k_grid, kp_grid]

# println("min = ", minimum(W_final_gal), ", max = ", maximum(W_final_gal))
# println("mean = ", mean(W_final_gal), ", std = ", std(W_final_gal))
# if  any(isnan, W_final_gal)   
#     println("There are NaN values in W_final_gal")
# end
# if any(isinf, W_final_gal)   
#     println("There are infinite values in W_final_gal")
# end
# if any(iszero, W_final_gal)
#     println("There are zero values in W_final_gal")
# end

println(minimum(k_grid), " to ", maximum(k_grid))
println(k_grid[1], " to ", k_grid[end])

idx = sortperm(k_grid)
idx_p = sortperm(kp_grid)

# genera tick automatici alle potenze di 10 nel range dei dati
xticks_vals = 10.0 .^ (floor(log10(minimum(k_grid))):ceil(log10(maximum(k_grid))))
yticks_vals = 10.0 .^ (floor(log10(minimum(kp_grid))):ceil(log10(maximum(kp_grid))))

heatmap(k_grid[idx], kp_grid[idx_p],
        W_final_gal[1,idx,idx_p]/maximum(W_final_gal[1,idx,idx_p]),
        title=L"W_{final}^{gg}"*"(at fixed "*L"ℓ)",
        xscale = :log10, yscale = :log10,
        xlabel=L"k", ylabel=L"k_p",
        c = c, xticks = xticks_vals, yticks = yticks_vals,
        yguidefontsize = yguidefontsize, xguidefontsize = xguidefontsize, titlefontsize = titlefontsize,
        dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin, size = size_heatmap
        )


savefig(joinpath(plot_subdir, "Sl_plots/k_grid_sorted_k_grid_sorted.png"))

println(size(W_tilde))

heatmap(1:Nk, 1:Nkp, W_final_gal[10,:,:]/maximum(W_final_gal[10,:,:]), 
    title = L"W_{final}^{gg}"*"(at fixed "*L"ℓ)", 
    xlabel = L"N_k", ylabel = L"N_{kp}", 
    colorbar = true,
    c = c, yguidefontsize = yguidefontsize, xguidefontsize = xguidefontsize, titlefontsize = titlefontsize,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin, size = size_heatmap
    )


savefig(joinpath(plot_subdir, "Sl_plots/Nk_Nkp.png"))

idx_p = sortperm(kp_grid)
yticks_vals = 10.0 .^ (floor(log10(minimum(k_grid))):ceil(log10(maximum(k_grid))))
heatmap(1:Nk, 
        kp_grid[idx_p], 
        W_final_gal[90, :, idx_p] / maximum(W_final_gal[90, :, idx_p]), 
        title = L"W_{final}^{gg}"*L"("*"at fixed "*L"ℓ = 90)", 
        xlabel = L"N_k", ylabel = L"\log_{10}(k_{grid})",
        yscale = :log10, colorbar = true,
        c = c, yticks = yticks_vals, yguidefontsize = yguidefontsize, xguidefontsize = xguidefontsize, titlefontsize = titlefontsize,
        dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin, size = size_heatmap
        )


savefig(joinpath(plot_subdir, "Sl_plots/Nk_kp_grid_sorted.png"))

#3D matter power spectrum
pk_dict = npzread("blast_code/data/pk.npz")
Pklin = pk_dict["pk_lin"]
Pknonlin = pk_dict["pk_nl"]
k_pk = pk_dict["k"]
z_pk = pk_dict["z"]
#Interpolating the power spectrum: Linear P(k) - Non-linear P(k)
y_pk = LinRange(log10(first(k_pk)),log10(last(k_pk)), length(k_pk))
x_pk = LinRange(first(z_pk), last(z_pk), length(z_pk))
InterpPmm = Interpolations.interpolate(log10.(Pklin),BSpline(Cubic(Line(OnGrid()))))
InterpPmm = scale(InterpPmm, (x_pk, y_pk))
InterpPmm = Interpolations.extrapolate(InterpPmm, Line())
InterpPmm_nl = Interpolations.interpolate(log10.(Pknonlin),BSpline(Cubic(Line(OnGrid()))))
InterpPmm_nl = scale(InterpPmm_nl, x_pk, y_pk)
InterpPmm_nl = Interpolations.extrapolate(InterpPmm_nl, Line())
power_spectrum(k_pk, χ1, χ2) = @. sqrt(10^InterpPmm(z_of_χ(χ1),log10(k_pk)) * 10^InterpPmm(z_of_χ(χ2),log10(k_pk)))
power_spectrum_nl(k_pk, χ1, χ2) = @. sqrt(10^InterpPmm_nl(z_of_χ(χ1),log10(k_pk)) * 10^InterpPmm_nl(z_of_χ(χ2),log10(k_pk)))

plot(k_pk, power_spectrum.(k_pk, 1000.0, 1000.0), 
     label="Linear P(k)", 
     xscale=:log10, 
     yscale=:log10, 
     title=L"$P(k)$ at $z=0$", titlefontsize=20,
     xlabel=L"$k \; (h/\mathrm{Mpc})$", 
     ylabel=L"$P(k) \; ((\mathrm{Mpc}/h)^3)$", 
     labelfontsize=15)
plot!(k_pk, power_spectrum_nl.(k_pk, 1000.0, 1000.0), 
      label="Non-linear P(k)", 
      xscale=:log10, 
      yscale=:log10,
      size=size_plot, legendfontsize=legendfontsize, dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
      )

Pk_grid = power_spectrum.(k_grid, kmin, kmax)
w_k = Blast.get_clencurt_weights(kmin, kmax, Nk)
weight_gal = w_k .* k_grid.^2 .* Pk_grid 
;

###alternative: I'd need to change also the integration of W_tilde
# log_k_grid = range(log(kmin), log(kmax), length=Nk)
# k_grid_2 = exp.(log_k_grid)
# dlnk = step(log_k_grid)
# w_log = ones(Nk) * dlnk
# w_log[1] /= 2.0
# w_log[end] /= 2.0
# weight_gal = w_log .* k_grid_2.^3 .* Pk_grid

abstract type AbstractProbe end
struct Galaxy <: AbstractProbe end
struct Shear <: AbstractProbe end
factorial_frac(ℓ) = (ℓ + 2.0) * (ℓ + 1.0) * ℓ * (ℓ - 1.0)
get_ell_prefactor(::Galaxy, ::Galaxy, ℓ) = @. (2 / π) * ones(length(ℓ))
get_ell_prefactor(::Galaxy, ::Shear,  ℓ) = @. (2 / π) * sqrt(factorial_frac(ℓ))
get_ell_prefactor(::Shear,  ::Galaxy, ℓ) = @. (2 / π) * sqrt(factorial_frac(ℓ))
get_ell_prefactor(::Shear,  ::Shear,  ℓ) = @. (2 / π) * factorial_frac(ℓ)
pref_gg = get_ell_prefactor(Galaxy(), Galaxy(), ℓ)
pref_gg = reduce(vcat, pref_gg)
pref_gs = get_ell_prefactor(Galaxy(), Shear(), ℓ)
pref_gs = reduce(vcat, pref_gs)
pref_gg = reduce(vcat, pref_gg)
pref_ss = get_ell_prefactor(Shear(), Shear(), ℓ)
pref_ss = reduce(vcat, pref_ss);

println("SIZES")
println("weight_gal -> ", size(weight_gal))
println("W_final_gal -> ", size(W_final_gal))
println("pref_gg -> ", size(pref_gg), "\npref_ss -> ", size(pref_ss), "\npref_gs -> ", size(pref_gs))

S_lkk_gg = zeros(Float64, size(W_final_gal, 3), size(W_final_gal, 3), length(ℓ))
@tullio S_lkk_gg[kp, kpp, li] = pref_gg[li] * weight_gal[k] * W_final_gal[li, k, kp] * W_final_gal[li, k, kpp]
npzwrite(joinpath(quantity_subdir, "Sl/S_lkk_gg.npy"), S_lkk_gg)
println("Size of S_l (kp, kpp) (gal-gal): \n(Nk, Nkp, NL) -> ", size(S_lkk_gg))
;

npzwrite(joinpath(quantity_subdir, "Sl/S_lkk_gg.npy"), S_lkk_gg)

xref = (xmax - xmin ) * 0.5
ℓ_to_k = ℓ -> ℓ ./ xref
ℓ_ticks = 1:50:200
k_ticks = ℓ_to_k.(ℓ_ticks)
k_ticklabels = [string(round(k, digits=3)) for k in k_ticks];

plot(ℓ, S_lkk_gg[1,1,:],
      color = colors[1],
      label = L"i, j=1")

plot!(xaxis = L"\ell",
      ylabel = L"S_\ell",
      xscale = :log10,
     )
plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xscale = :log10,
      xlims = (minimum(ℓ), maximum(ℓ)),
      label = ""
     )
plot!(label = "Beyond BLAST", 
      size = size_plot,
      title = L"S_\ell^{gg} (i = j)",
      titlefontsize = 20,
      titleposition = :left,
      dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin, rightmargin = 5Plots.mm, topmargin = 5Plots.mm)

plot(k_grid, S_lkk_gg[:,1,1],
      color = colors[1],
      label = L"i, j=1")

plot!(xaxis = L"k_{grid} \; (h/\mathrm{Mpc})",
      ylabel = L"S_\ell",
      #xscale = :log10,
     )
plot!(label = "Beyond BLAST", 
      size = size_plot,
      title = L"S_\ell^{gg} (k, kp = 1, \ell = 20)",
      titlefontsize = 20,
      titleposition = :left,
      dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin, rightmargin = 5Plots.mm, topmargin = 5Plots.mm)

plt = plot()

for i in 1:Nkp
    plot!(plt, ℓ, S_lkk_gg[i,i,:], line_z = i,
          label = L"i = %$(i)",
          color = colors, linewidth = 2)
end

plot!(plt, label="Beyond BLAST", 
     xaxis=L"\ell", ylabel=L"S_\ell", 
     xscale = :log10
    )
plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ)),
      xscale = :log10,
      label = ""
     )
plot!(legend = false, colorbar = true, colorbar_title = L"i = j",
      clims = (1, Nkp),
      titleposition = :left,
      title=L"S_\ell^{gg} (i = j)",
      leftmargin = leftmargin, size=size_plot, titlefontsize=titlefontsize, bottommargin = bottommargin, dpi = dpi
     )
plt

plt = plot()

for i in 1:Nk
    plot!(plt, ℓ, S_lkk_gg[i,i,:]/maximum(S_lkk_gg[i,i,:]), line_z = i,
          label = L"i = %$(i)",
          color = colors, linewidth = 1)
end

plot!(plt, 
    clims = (1, Nk),
    xaxis=L"\ell", yaxis=L"S_\ell", 
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ)),
      label = ""
     )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j",
    title=L"S_\ell^{gg} (i = j)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt

savefig(joinpath(plot_subdir, "Sl_plots/S_lkk_gg_fixed_l_various.png"))

plt = plot()

a = 1:5:Nk
for i in a
    plot!(plt, ℓ, S_lkk_gg[i,i,:]/maximum(S_lkk_gg[i,i,:]), line_z = i,
          label = L"i = %$(i)", color = colors, linewidth = 1
          )
end

plot!(plt, 
    xaxis=L"\ell", yaxis=L"S_\ell", 
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ))     
      )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"i_k", clims = (1, Nkp),
    title=L"S_\ell^{gg} (i = j)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt

i_fixed = 1
for j in 1:Nkp
    plot!(plt, ℓ, S_lkk_gg[i_fixed, j, :],
          line_z = j, color = colors, linewidth = 2
        )                  
end

plot!(plt,
    xaxis = L"\ell", yaxis = L"S_\ell",
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ))
     )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j", clims = (1, Nkp),
    title = L"S_\ell^{gg} (i = 1,\ j)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt


plt = plot()
jj = 96
i_fixed = 1
for j in 1:jj
    plot!(plt, ℓ, S_lkk_gg[i_fixed, j, :] / maximum(S_lkk_gg[i_fixed, j, :]),
          line_z = j, color = colors, linewidth = 2
          )
end

plot!(plt,
    xaxis = L"\ell", yaxis = L"S_\ell",
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ))
     )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j", clims = (1, Nkp),
    title = L"S_\ell^{gg} (i = 1,\ j)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt

plt = plot()
jj = 100
i_fixed = 50
for j in 1:jj
    plot!(plt, ℓ, S_lkk_gg[i_fixed, j, :],
          line_z = j, color = colors,
          linewidth = 2
          )
end

plot!(plt,
    xaxis = L"\ell", yaxis = L"S_\ell",
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ))
    )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j", clims = (1, Nkp),
    title = L"S_\ell^{gg} (i = 10,\ j)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt

plot(
    ℓ,
    S_lkk_gg[1, 1, :] .* ℓ .* (ℓ .+ 1),
    xaxis = L"\ell",
    yaxis = L"\ell(\ell+1)S_\ell",
)

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ)),
      label = ""
     )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j", clims = (1, Nkp),
    title = L"\ell(\ell+1)S_\ell^{gg}(k_1,k_2)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt = plot()

for j in 1:Nkp
    plot!(plt, ℓ, S_lkk_gg[j, j, :] .* ℓ .* (ℓ .+ 1),
          line_z = j, color = colors, linewidth = 2
          )
end

plot!(plt,
    xaxis = L"\ell", 
    yaxis = L"\ell(\ell+1)S_\ell",
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ))
     )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j", clims = (1, Nkp),
    title = L"\ell(\ell+1)S_\ell^{gg}(k_1,k_2)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt

plt = plot()

for j in a
    plot!(plt, ℓ, S_lkk_gg[j, j, :] .* ℓ .* (ℓ .+ 1),
          line_z = j, color = colors, linewidth = 2
          )                  
end

plot!(plt,
    xaxis = L"\ell", 
    yaxis = L"\ell(\ell+1)S_\ell",
    xscale = :log10
    )

plot!(twiny(),
      xaxis = L"k \; (h/\mathrm{Mpc})",
      xticks = (ℓ_ticks, k_ticklabels),
      xlims = (minimum(ℓ), maximum(ℓ))
     )

plot!(label = "Beyond BLAST", 
    legend = false, colorbar = true, colorbar_title = L"j", clims = (1, Nkp),
    title = L"\ell(\ell+1)S_\ell^{gg}(k_1,k_2)",
    size=size_plot, titlefontsize=20, 
    titleposition = :left,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
    )

plt

println(size(S_lkk_gg))

diagS = [diag(S_lkk_gg[:, :, i]) for i in 1:size(S_lkk_gg, 3)]

println("Size of S_lkk_gg: ", size(S_lkk_gg))
println("Size of diagS: ", size(diagS))
println("Size of k_grid: ", size(k_grid))
println("Size of kp_grid: ", size(kp_grid))

nshow = 100
idx = round.(Int, range(1, size(S_lkk_gg, 3), length=nshow))

p = plot(
    xlabel = L"k \; (h/\mathrm{Mpc})",
    ylabel = L"S_\ell(k,k)",
    title  = L"S_\ell^{gg}(k,k)\ \mathrm{for\ different}\ \ell\ \mathrm{values}",
    colorbar_title = L"j",
    clims = (1, Nkp),                   
    legend = false,
    colorbar = true,
    lw = 2,
    size = size_plot, titlefontsize = titlefontsize,
    dpi = dpi, leftmargin = leftmargin, bottommargin = bottommargin
)

for ii in idx
    plot!(p, kp_grid, diagS[ii]/maximum(diagS[ii]), line_z = ii, color = colors;
        label = L"\ell = %$(ℓ[ii])")
end

p
