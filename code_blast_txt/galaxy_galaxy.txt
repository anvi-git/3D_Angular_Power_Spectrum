module galaxy_galaxy

using LaTeXStrings
using Plots
using Cosmology
using QuadGK
using DataInterpolations
using NPZ

import Cosmology: AbstractCosmology
import Main.Blast
import PhysicalConstants.CODATA2018: c_0 # speed of light in m/s
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

function heath_integral(cosmo, z)
    integrand(zp) = (1.0 + zp) / (Blast.compute_adimensional_hubble_factor(zp, cosmo)^3)
    integral_val, _ = quadgk(integrand, z, Inf, rtol=1e-10)
    Ez = Blast.compute_adimensional_hubble_factor(z, cosmo)
    return Ez * integral_val
end
function compute_growth_factor(cosmo, z)
    D_unnorm_z0 = heath_integral(cosmo, 0.0)
    D_array = [heath_integral(cosmo, zz) / D_unnorm_z0 for zz in z]
    return D_array
end

function galaxy_prefactor(
    x::AbstractVector,
    z::AbstractVector,
    cosmo;
    output_dir::AbstractString
)
    # setting the output directory for plots
    println("Plots will be saved in: ", joinpath(output_dir, "plots/kernels/"))
    println("Quantities will be saved in: ", joinpath(output_dir, "quantities/"))
    # computing the plot comoving distance vs redsfhit
    plot(z, x, label = L"z(\chi)", xlabel = L"z", 
         ylabel = L"\chi [Mpc/h]", legend = :topleft,
         size=(800,600), dpi = 200)
    savefig(joinpath(output_dir, "plots/kernels/chi_vs_z.png"))
    println("Comoving distance vs redshift plot saved.")
    ######
    # BIAS 
    # compute the bias. the equation is b(z) = b_0 * sqrt(1+z). we set b_0 = 1.0.!
    b_0 = 1.0
    # bz_array = b_0 .* sqrt.(1 .+ z); # bias as a function of redshift
    bz_array = @. 0.98 + 1.24 * z - 1.72 * z^2 + 1.28 * z^3    
    #plot
    plot(z, bz_array, label = "bias", xlabel = L"z", 
         ylabel = L"b(z)", legend = :topleft,
         size=(800,600), dpi = 200)
    savefig(joinpath(output_dir, "plots/kernels/bias.png"))
    println("Bias plot saved.")
    ######
    # GROWTH FACTOR
    D_growth_array = compute_growth_factor(cosmo, z)
    #plot
    plot(z, D_growth_array, label = "growth array", xlabel = L"z", 
         ylabel = L"D(z)", legend = :topleft,
         size=(800,600), dpi = 200)
    savefig(joinpath(output_dir, "plots/kernels/growth.png"))
    println("Growth factor plot saved.")
    ######
    # Hubble factor
    Hubble_array = Blast.compute_hubble_factor.(z, Ref(cosmo)) ./ C_LIGHT
    #plot
    plot(z, Hubble_array, label = "Hubble factor", xlabel = L"z", 
    ylabel = L"H(z) [Mpc/h]/c", legend = :topleft,
         size=(800,600), dpi = 200)
    savefig(joinpath(output_dir, "plots/kernels/hubble.png"))
    println("Hubble factor plot saved.")
    ######
    #n(z)  
    z_0 = 0.9/sqrt(2)
    A = 1.5/z_0
    alpha = 2
    beta = 1.5
    z_of_χ = DataInterpolations.AkimaInterpolation(z, x); # z(χ)
    dz_dchi = [DataInterpolations.derivative(z_of_χ, chi) for chi in x]
    nz = A .* (z / z_0).^alpha .* exp.(-(z / z_0).^beta) # redshift distribution of sources, normalized to 1
    nz_interp = DataInterpolations.AkimaInterpolation(nz, z, extrapolation=ExtrapolationType.Linear)
    integral_nz, _ = quadgk(x -> nz_interp(x), minimum(z), maximum(z), rtol=1e-10)
    nz_norm = nz ./ integral_nz
    n_chi_norm = nz_norm .* dz_dchi    
    ### CASO n(z) con picco su z = 1
    # z_0 = 1.0
    # alpha = 10.0
    # beta  = 10.0
    # nz = (z ./ z_0).^alpha .* exp.(- (z ./ z_0).^beta)
    # nz_norm = nz ./ quadgk(x -> DataInterpolations.AkimaInterpolation(nz, z, extrapolation=ExtrapolationType.Linear)(x),
    #                     minimum(z), maximum(z))[1]
    #plot
    plot(z, nz, label=L"n(z)", xlabel=L"z", 
         ylabel=L"n(z)")
    savefig(joinpath(output_dir, "plots/kernels/nz.png"))
    println("n(z) vs redshift plot saved.")
    plot(z, n_chi_norm, label=L"n(z) normalized", xlabel=L"z", 
          ylabel=L"n(z) normalized")
    savefig(joinpath(output_dir, "plots/kernels/nz_normalized.png"))
    println("Normalized n(z) vs redshift plot saved.")
    plot(z, nz, label=L"n(z)", xlabel=L"z", 
         ylabel=L"n(z)")
    plot!(z, n_chi_norm, label=L"n(z) normalized", xlabel=L"z", 
          ylabel=L"n(z) normalized")
    savefig(joinpath(output_dir, "plots/kernels/nz_both.png"))
    println("n(z) vs redshift plot saved.")
    ######
    # W(z)
    prefac = x.^2 .* bz_array .* D_growth_array .* n_chi_norm .* Hubble_array
    #plot
    plot(z, 
          prefac, 
          label=L"W(z) = \frac{n(z) b(z) D(z) \chi(z)^2 H(z)}{c}", 
          xlabel=L"z", 
          ylabel=L"W(z) [Mpc/h]^2", 
          title=L"W(z)", 
          legend=:topright, size=(800,600), titlefontsize=15, labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/kernels/Wz.png"))
    println("W(z) plot saved.")
    #normalized plot
    plot(z, 
          bz_array ./ maximum(bz_array), 
          label=L"b(z)", 
          xlabel=L"z", 
          ylabel=L"Normalized b(z)", 
          title=L"Normalized b(z)",
                size=(800,600), dpi = 200)
        plot!(z, 
                x ./ maximum(x), 
                label=L"\chi(z)",
                title = L"Normalized \chi",
                xlabel=L"z", 
                ylabel=L"Normalized \chi",
                size=(800,600), dpi = 200)
        plot!(z, 
                D_growth_array ./ maximum(D_growth_array), 
                label=L"D(z)", 
                xlabel=L"z", 
                ylabel=L"Normalized D(z)", 
                title=L"Normalized D(z)",
                legend=:bottomright, size=(800,600), titlefontsize=15, labelfontsize=15, legendfontsize=15, dpi = 200)
        plot!(z, 
                n_chi_norm ./ maximum(n_chi_norm), 
                label=L"n(z)", 
                xlabel=L"z", 
                ylabel=L"Normalized n(z)", 
                title=L"Normalized n(z)",
                size=(800,600), dpi = 200)
        plot!(z, 
                  Hubble_array ./ maximum(Hubble_array), 
                  label=L"H(z)", 
                  xlabel=L"z", 
                  ylabel=L"Normalized H(z)", 
                  title=L"Normalized H(z)", 
                  size=(800,600), dpi = 200)
        plot!(z, 
                prefac ./ maximum(prefac), 
                label=L"\frac{n(z) b(z) D(z) \chi(z)^2 H(z)}{c}", 
                xlabel=L"z", 
                ylabel="Normalized Kernel", 
                title="Normalized Kernel", 
                legend=:bottomright, size=(800,600), titlefontsize=15, labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/kernels/normalized_kernel_galaxy.png"))
    println("Normalized kernel plot saved.")

    plot(z, 
        prefac, 
        label=L"\frac{n(z) b(z) D(z) \chi(z)^2 H(z)}{c}", 
        xlabel=L"z", 
        ylabel="Unnormalized Kernel", 
        title="Unnormalized Kernel", 
        legend=:bottomright, size=(800,600), titlefontsize=15, labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/kernels/unnormalized_kernel_galaxy.png"))
    println("Unnormalized kernel plot saved.")

    npzwrite(joinpath(output_dir, "quantities/prefac_galaxy.npy"), prefac)
    println("W(z) prefactor saved: ", joinpath(output_dir, "quantities/prefac_galaxy.npy"))
    return prefac, bz_array, D_growth_array, Hubble_array, n_chi_norm
end

function galaxy_prefactor_cheb(xmin::Number, 
                               xmax::Number, 
                               n_cheb::Int,
                               z::AbstractVector, 
                               x::AbstractVector, 
                               bias::AbstractVector, 
                               growth::AbstractVector, 
                               Hubble_param::AbstractVector, 
                               n_chi_norm::AbstractVector; 
                               output_dir::AbstractString)

    chi_cheb_nodes = Blast.get_clencurt_grid(xmin, xmax, n_cheb)

    z_interp = DataInterpolations.AkimaInterpolation(z, x, extrapolation=ExtrapolationType.Linear)
    z_cheb_nodes = z_interp.(chi_cheb_nodes) 

    b_interp = DataInterpolations.AkimaInterpolation(bias, x, extrapolation=ExtrapolationType.Linear)
    b_cheb = b_interp.(chi_cheb_nodes)

    D_interp = DataInterpolations.AkimaInterpolation(growth, x, extrapolation=ExtrapolationType.Linear)
    D_cheb = D_interp.(chi_cheb_nodes)

    nz_interps = DataInterpolations.AkimaInterpolation(n_chi_norm, x, extrapolation=ExtrapolationType.Linear)
    nz_cheb = nz_interps.(chi_cheb_nodes)

    H_interp = DataInterpolations.AkimaInterpolation(Hubble_param, x, extrapolation=ExtrapolationType.Linear)
    Hubble_cheb = H_interp.(chi_cheb_nodes)

    W_vals = zeros(length(chi_cheb_nodes))
    W_vals = chi_cheb_nodes.^2 .* b_cheb .* D_cheb .* nz_cheb .* Hubble_cheb

    npzwrite(joinpath(output_dir, "quantities/prefac_on_cheb_galaxy.npy"), W_vals)
    println("W(z) prefactor on Chebyshev grid saved: ", joinpath(output_dir, "quantities/prefac_on_cheb_galaxy.npy"))
    return W_vals
end

function compute_prefactor_chebcoeffs(W_vals::AbstractVector; output_dir::AbstractString)
    plan = Blast.plan_fft(W_vals)
    cheb_coeff = Blast.fast_chebcoefs(W_vals, plan)
    npzwrite(joinpath(output_dir, "quantities/W_cheb_coeff_galaxy.npy"), cheb_coeff)
    println("W(z) Chebyshev coefficients saved: ", joinpath(output_dir, "quantities/W_cheb_coeff_galaxy.npy"))

    return cheb_coeff

end 
end # module



