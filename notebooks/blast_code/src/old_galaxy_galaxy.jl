module galaxy_galaxy

using LaTeXStrings
using Plots
using Cosmology
using QuadGK
using DataInterpolations
using NPZ

import Cosmology: AbstractCosmology
import Main.Blast
import PhysicalConstants.CODATA2018: c_0
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

function heath_integral(cosmo, z)
    integrand(zp) = (1.0 + zp) / (Blast.compute_adimensional_hubble_factor(zp, cosmo)^3)
    integral_val, _ = quadgk(integrand, z, Inf, rtol=1e-10)
    return Blast.compute_adimensional_hubble_factor(z, cosmo) * integral_val
end
function compute_growth_factor(cosmo, z_range)
    D_unnorm_z0 = heath_integral(cosmo, 0.0)
    D_array = [heath_integral(cosmo, zz) / D_unnorm_z0 for zz in z_range]
    return D_array
end

function galaxy_prefactor(
    x_range::AbstractVector,
    z_range::AbstractVector,
    cosmo;
    output_dir::AbstractString
)
    # computing the plot comoving distance vs redsfhit
    plot(z_range, x_range, label = "z(χ)", xlabel = L"z", 
         ylabel = L"\chi [Mpc/h]", legend = :topleft,
         size=(800,600))
    savefig(joinpath(output_dir, "plots/kernels/chi_vs_z.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/chi_vs_z.png"))
    ######
    # BIAS 
    # compute the bias. the equation is b(z) = b_0 * sqrt(1+z). we set b_0 = 1.0.!
    b_0 = 1.0
    bz_array = zeros(length(z_range))
    bz_array = b_0 .* sqrt.(1 .+ z_range); # bias as a function of redshift
    #plot
    plot(z_range, bz_array, label = "bias", xlabel = L"z", 
         ylabel = L"b(z)", legend = :topleft,
         size=(800,600))
    savefig(joinpath(output_dir, "plots/kernels/bias.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/bias.png"))
    ######
    # GROWTH FACTOR
    D_growth_array = compute_growth_factor(cosmo, z_range)
    #plot
    plot(z_range, D_growth_array, label = "growth array", xlabel = L"z", 
         ylabel = L"D(z)", legend = :topleft,
         size=(800,600))
    savefig(joinpath(output_dir, "plots/kernels/growth.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/growth.png"))
    ######
    # Inverse Hubble factor
    inv_Hubble_array = C_LIGHT ./ Blast.compute_hubble_factor.(z_range, Ref(cosmo))
    #plot
    plot(z_range, inv_Hubble_array, label = "inverse Hubble factor", xlabel = L"z", 
    ylabel = L"c/H(z) [Mpc/h]", legend = :topleft,
         size=(800,600))
    savefig(joinpath(output_dir, "plots/kernels/inv_hubble.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/inv_hubble.png"))
    ######
    #n(z)  
    z_0 = 0.9/sqrt(2)
    A = 1.5/z_0
    alpha = 2
    beta = 1.5
    nz = A .* (z_range / z_0).^alpha .* exp.(-(z_range / z_0).^beta) # redshift distribution of sources, normalized to 1
    nz_norm = nz ./ quadgk(x -> DataInterpolations.AkimaInterpolation(nz, z_range, extrapolation=ExtrapolationType.Linear)(x),
                      minimum(z_range), maximum(z_range))[1]
    ### CASO n(z) con picco su z = 1
    # z_0 = 1.0
    # alpha = 10.0
    # beta  = 10.0
    # nz = (z_range ./ z_0).^alpha .* exp.(- (z_range ./ z_0).^beta)
    # nz_norm = nz ./ quadgk(x -> DataInterpolations.AkimaInterpolation(nz, z_range, extrapolation=ExtrapolationType.Linear)(x),
    #                     minimum(z_range), maximum(z_range))[1]
    #plot
    plot(z_range, nz, label="n(z)", xlabel="z", 
         ylabel="n(z)", title="Redshift distribution of sources",
         size=(800,600))
    plot!(z_range, nz_norm, label="n(z) normalized", xlabel="z", 
          ylabel="n(z)", title="Redshift distribution of sources",
          size=(800,600))
    savefig(joinpath(output_dir, "plots/kernels/nz.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/nz.png"))
    ######
    # W(z)
    prefac = x_range.^2 .* bz_array .* D_growth_array .* nz_norm
    #plot
    plot(z_range, 
          prefac, 
          label=L"$ W(z) = n(z) b(z) D(z) \chi(z)^2$", 
          xlabel=L"$z$", 
          ylabel=L"$W(z)$ $[Mpc/h]^2$", 
          title=L"$W(z)$", 
          legend=:topright, size=(800,600), titlefontsize=15, labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/kernels/Wz.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/Wz.png"))
    #normalized plot
    plot(z_range, 
          bz_array ./ maximum(bz_array), 
          label=L"$b(z)$", 
          xlabel=L"$z$", 
          ylabel=L"Normalized $b(z)$", 
          title=L"Normalized $b(z)$",
                size=(800,600))
        plot!(z_range, 
                x_range ./ maximum(x_range), 
                label=L"$\chi$",
                title = L"Normalized $\chi$",
                xlabel=L"$z$", 
                ylabel=L"Normalized $\chi$",
                size=(800,600))
        plot!(z_range, 
                D_growth_array ./ maximum(D_growth_array), 
                label=L"$D(z)$", 
                xlabel=L"$z$", 
                ylabel="Normalized quantities", 
                title="Normalized quantities",
                legend=:bottomright, size=(800,600), titlefontsize=15, labelfontsize=15, legendfontsize=15, dpi = 200)
        plot!(z_range, 
                nz_norm ./ maximum(nz_norm), 
                label=L"$n(z)$", 
                xlabel=L"$z$", 
                ylabel=L"Normalized $n(z)$", 
                title=L"Normalized $n(z)$",
                size=(800,600))
        # plot!(z_range, 
        #           inv_Hubble_array ./ maximum(inv_Hubble_array), 
        #           label=L"$1/H(z)$", 
        #           xlabel=L"$z$", 
        #           ylabel=L"Normalized $1/H(z)$", 
        #           title=L"Normalized $1/H(z)$")
        plot!(z_range, 
                prefac ./ maximum(prefac), 
                label=L"$n(z) b(z) D(z) \chi(z)^2$", 
                xlabel=L"$z$", 
                ylabel="Normalized Kernel", 
                title="Normalized Kernel", 
                legend=:bottomright, size=(1400,800), titlefontsize=15, labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/kernels/normalized_kernel_galaxy.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/normalized_kernel_galaxy.png"))
    plot(z_range, 
        prefac, 
        label=L"$n(z) b(z) D(z) \chi(z)^2$", 
        xlabel=L"$z$", 
        ylabel="Unnormalized Kernel", 
        title="Unnormalized Kernel", 
        legend=:bottomright, size=(1400,800), titlefontsize=15, labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/kernels/unnormalized_kernel_galaxy.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/kernels/unnormalized_kernel_galaxy.png"))
    npzwrite(joinpath(output_dir, "quantities/prefac_galaxy.npy"), prefac)
    println("Wrote W(z) prefactor to: ", joinpath(output_dir, "quantities/prefac_galaxy.npy"))
    return prefac, bz_array, D_growth_array, inv_Hubble_array, nz_norm
end

function galaxy_prefactor_cheb(zmin::Number, 
                               zmax::Number, 
                               n_cheb::Int,
                               z_range::AbstractVector, 
                               chi::AbstractVector, 
                               bias::AbstractVector, 
                               growth::AbstractVector, 
                               inv_Hubble::AbstractVector, 
                               nz_norm::AbstractVector; 
                               output_dir::AbstractString)

    z_cheb_nodes = Blast.get_clencurt_grid_z(zmin, zmax, n_cheb)
    chi_interp = DataInterpolations.AkimaInterpolation(chi, z_range, extrapolation=ExtrapolationType.Linear)
    chi_cheb = chi_interp.(z_cheb_nodes)
    b_interp = DataInterpolations.AkimaInterpolation(bias, z_range, extrapolation=ExtrapolationType.Linear)
    b_cheb = b_interp.(z_cheb_nodes)
    D_interp = DataInterpolations.AkimaInterpolation(growth, z_range, extrapolation=ExtrapolationType.Linear)
    D_cheb = D_interp.(z_cheb_nodes)
    nz_interps = DataInterpolations.AkimaInterpolation(nz_norm, z_range, extrapolation=ExtrapolationType.Linear)
    nz_cheb = nz_interps.(z_cheb_nodes)
    W_vals = zeros(length(z_cheb_nodes))
    W_vals .= chi_cheb.^2 .* b_cheb .* D_cheb .* nz_cheb
    npzwrite(joinpath(output_dir, "quantities/prefac_on_cheb_galaxy.npy"), W_vals)
    println("Wrote W(z) prefactor on Chebyshev grid to: ", joinpath(output_dir, "quantities/prefac_on_cheb_galaxy.npy"))
    return W_vals
end

function compute_prefactor_chebcoeffs(W_vals::AbstractVector; output_dir::AbstractString)
    plan = Blast.plan_fft(W_vals)
    cheb_coeff = Blast.fast_chebcoefs(W_vals, plan)
    npzwrite(joinpath(output_dir, "quantities/W_cheb_coeff_galaxy.npy"), cheb_coeff)
    println("Wrote W(z) Chebyshev coefficients to: ", joinpath(output_dir, "quantities/W_cheb_coeff_galaxy.npy"))

    return cheb_coeff

end 
end # module



