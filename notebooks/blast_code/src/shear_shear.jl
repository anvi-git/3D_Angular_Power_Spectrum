module shear_shear

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

function shear_prefactor(
    x_range::AbstractVector,
    z_range::AbstractVector,
    cosmo;
    output_dir::AbstractString
)
    # computing comoving distance 
    com_dist = Blast.compute_χ.(z_range, Ref(cosmo));
    plot(z_range, com_dist, label = "z(χ)", xlabel = L"z", ylabel = L"\chi [Mpc/h]", legend = :topleft)
    savefig(joinpath(output_dir, "plots/chi_vs_z_shear.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/chi_vs_z_shear.png"))
    ######
    # n(z)  
    z_0 = 0.9/sqrt(2)
    alpha = 2
    beta = 1.5
    nz = (z_range / z_0).^alpha .* exp.(-(z_range / z_0).^beta) # redshift distribution of sources, normalized to 1
    nz_norm = nz ./ quadgk(x -> DataInterpolations.AkimaInterpolation(nz, z_range, extrapolation=ExtrapolationType.Linear)(x),
                      minimum(z_range), maximum(z_range))[1]
    nz_interp = DataInterpolations.AkimaInterpolation(nz_norm, z_range, extrapolation=ExtrapolationType.Linear)
    #plot
    plot(z_range, nz, label="n(z)", xlabel="z", ylabel="n(z)", title="Redshift distribution of sources")
    plot!(z_range, nz_norm, label="n(z) normalized", xlabel="z", ylabel="n(z)", title="Redshift distribution of sources")
    savefig(joinpath(output_dir, "plots/nz_shear.png"))
    println("Plot saved: ", joinpath(output_dir, "plots/nz_shear.png"))
    prefac_shear = similar(z_range, Float64)
    pref = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    zmax_int = maximum(z_range)
    for i in eachindex(z_range)
        zi = z_range[i]
        χi = com_dist[i]
        integrand(zs) = begin
            χs = Blast.compute_χ(zs, cosmo)
            χs <= χi ? 0.0 : nz_interp(zs) * (1.0 - χi / χs)
        end
        lens_int, _ = quadgk(integrand, zi, zmax_int)
        prefac_shear[i] = pref * χi * (1.0 + zi) * lens_int
    end
    plot(z_range, prefac_shear,
         label="W(z)", xlabel="z", ylabel="W(z)",
         title="Prefactor for shear-shear correlation",
         legend=:topright, size=(800, 600), titlefontsize=15,
         labelfontsize=15, legendfontsize=7, dpi=200)
    savefig(joinpath(output_dir, "plots/quantities_shear.png"))
    plot(z_range, prefac_shear ./ maximum(prefac_shear),
         label="W(z) normalized", xlabel="z", ylabel="W(z)/max(W(z))",
         title="Normalized shear prefactor")
    savefig(joinpath(output_dir, "plots/normalized_quantities_shear.png"))
    npzwrite(joinpath(output_dir, "quantities/prefac_shear.npy"), prefac_shear)
    return prefac_shear
end

function shear_prefactor_cheb(
    zmin::Number,
    zmax::Number,
    n_cheb::Int,
    z_range::AbstractVector,
    prefac_shear::AbstractVector;
    output_dir::AbstractString
)
    z_cheb_nodes = Blast.get_clencurt_grid_z(zmin, zmax, n_cheb)
    W_interp = DataInterpolations.AkimaInterpolation(
        prefac_shear, z_range,
        extrapolation=ExtrapolationType.Linear
    )

    W_vals = W_interp.(z_cheb_nodes)

    npzwrite(joinpath(output_dir, "quantities/prefac_on_cheb_shear.npy"), W_vals)
    println("Wrote W(z) shear prefactor on Chebyshev grid to: ",
            joinpath(output_dir, "quantities/prefac_on_cheb_shear.npy"))

    return W_vals
end

function compute_prefactor_chebcoeffs(W_vals::AbstractVector; output_dir::AbstractString)
    plan = Blast.plan_fft(W_vals)
    cheb_coeff = Blast.fast_chebcoefs(W_vals, plan)
    npzwrite(joinpath(output_dir, "quantities/W_cheb_coeff_shear.npy"), cheb_coeff)
    println("Wrote W(z) Chebyshev coefficients to: ", joinpath(output_dir, "quantities/W_cheb_coeff_shear.npy"))

    return cheb_coeff

end 
end # module



