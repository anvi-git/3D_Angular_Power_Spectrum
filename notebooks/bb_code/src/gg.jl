import Cosmology: AbstractCosmology
import PhysicalConstants.CODATA2018: c_0 # speed of light in m/s
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

function heath_integral(cosmo, z; output_dir::Union{AbstractString, Nothing} = nothing)
    integrand(zp) = (1.0 + zp) / (bb.compute_adimensional_hubble_factor(zp, cosmo)^3)
    integral_val, _ = quadgk(integrand, z, Inf, rtol=1e-10)
    Ez = bb.compute_adimensional_hubble_factor(z, cosmo)
    res = Ez * integral_val

    return res
end

function compute_growth_factor(cosmo, z; output_dir::Union{AbstractString, Nothing} = nothing)
    if !isnothing(output_dir)
        append_to_log(output_dir, "=== Starting growth factor computation ===")
        append_to_log(output_dir, "Evaluating unnormalized growth factor at z = 0...")
    end
    D_unnorm_z0 = heath_integral(cosmo, 0.0; output_dir = output_dir)
    D_array = [heath_integral(cosmo, zz; output_dir = output_dir) / D_unnorm_z0 for zz in z]
    return D_array
end

function compute_Wx(
    x::AbstractVector,
    z::AbstractVector,
    cosmo;
    output_dir::AbstractString,
    plot_style::AbstractDict = Dict()
    )
    # setting the output directory for plots
    bb.append_to_log(output_dir, "=== Galaxy prefactor function started ===")
    bb.append_to_log(output_dir, "Processing galaxy prefactor calculation...")
    # computing the plot comoving distance vs redshift
    bb.append_to_log(output_dir, "=== Comoving distance vs redshift ===")
    plot(z, x, label = L"z(\chi)", xlabel = L"z", ylabel = L"\chi [Mpc/h]", legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/χ.png"))
    bb.append_to_log(output_dir, "Size of the comoving distance array x is $(length(x)).")
    bb.append_to_log(output_dir, "Done. Moving on to the bias...")

    ######
    # BIAS 
    bb.append_to_log(output_dir, "=== Bias ===")
    bb.append_to_log(output_dir, "Computing the bias...")
    b_0 = 1.0
    bz_array = @. 0.98 + 1.24 * z - 1.72 * z^2 + 1.28 * z^3    

    plot(z, bz_array, label = "bias", xlabel = L"z", ylabel = L"b(z)", legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/bias.png"))
    bb.append_to_log(output_dir, "Size of the bias array is $(length(bz_array)).")

    ######
    # GROWTH FACTOR
    bb.append_to_log(output_dir, "=== Growth factor ===")
    bb.append_to_log(output_dir, "Computing the growth factor...")
    D_growth_array = compute_growth_factor(cosmo, z; output_dir = output_dir)

    plot(z, D_growth_array, label = "growth array", xlabel = L"z", ylabel = L"D(z)", legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/growth.png"))
    bb.append_to_log(output_dir, "Size of the growth factor array is $(length(D_growth_array)).")

    ######
    # HUBBLE FACTOR
    bb.append_to_log(output_dir, "=== Hubble factor ===")
    bb.append_to_log(output_dir, "Computing the Hubble factor...")
    Hubble_array = bb.compute_hubble_factor.(z, Ref(cosmo)) ./ C_LIGHT

    plot(z, Hubble_array, label = "Hubble factor", xlabel = L"z", ylabel = L"H(z) [Mpc/h]/c", legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/h_over_c.png"))
    bb.append_to_log(output_dir, "Size of the Hubble factor array is $(length(Hubble_array)).")

    ######
    # NUMBER DENSITY n(z)
    bb.append_to_log(output_dir, "=== Number density ===")
    bb.append_to_log(output_dir, "Computing the number density...")
    z_0 = 0.9 / sqrt(2)
    A = 1.5 / z_0
    alpha = 2
    beta = 1.5
    z_of_χ = DataInterpolations.AkimaInterpolation(z, x)
    dz_dchi = [DataInterpolations.derivative(z_of_χ, chi) for chi in x]
    nz = A .* (z / z_0).^alpha .* exp.(-(z / z_0).^beta)
    nz_interp = DataInterpolations.AkimaInterpolation(nz, z, extrapolation=ExtrapolationType.Linear)
    integral_nz, _ = quadgk(x -> nz_interp(x), minimum(z), maximum(z), rtol=1e-10)
    nz_norm = nz ./ integral_nz
    n_chi_norm = nz_norm .* dz_dchi 

    plot(z, nz, label=L"n(z)", xlabel=L"z", ylabel=L"n(z)", legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/nz.png"))
    bb.append_to_log(output_dir, "Plotting unnormalized and normalized distributions...")

    plot(z, n_chi_norm, label=L"n(z) normalized", xlabel=L"z", ylabel=L"n(z) normalized", legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/nz_normalized.png"))
    bb.append_to_log(output_dir, "Normalized n(z) vs redshift plot saved.")

    plot(z, nz, label=L"n(z)", xlabel=L"z", ylabel=L"n(z)")
    plot!(z, n_chi_norm, label=L"n(z) normalized", xlabel=L"z", ylabel=L"n(z) normalized")
    plot!(legend = :topleft; plot_style...)
    savefig(joinpath(output_dir, "plots/W/nz_both.png"))
    bb.append_to_log(output_dir, "Combined n(z) plot saved.")

    bb.append_to_log(output_dir, "Size of the number density array is $(length(nz)).")
    bb.append_to_log(output_dir, "Size of the normalized number density array is $(length(n_chi_norm)).")

    ######
    # PREFACTOR W(z)
    bb.append_to_log(output_dir, "=== W(χ) ===")
    bb.append_to_log(output_dir, "Computing W(χ)...")
    W_x = x.^2 .* bz_array .* D_growth_array .* n_chi_norm .* Hubble_array

    plot(z, W_x, label=L"W(z) = \frac{n(z) b(z) D(z) \chi(z)^2 H(z)}{c}", xlabel=L"z", ylabel=L"W(z) [Mpc/h]^2", title=L"W(z)", legend=:topright; plot_style...)
    savefig(joinpath(output_dir, "plots/W/Wz.png"))
    bb.append_to_log(output_dir, "W(z) plot saved.")

    # Normalized comparison plot
    plot(z, bz_array ./ maximum(bz_array), label=L"b(z)", xlabel=L"z", ylabel=L"Normalized b(z)", title=L"Normalized b(z)"; plot_style...)
    plot!(z, x ./ maximum(x), label=L"\chi(z)", title = L"Normalized \chi", xlabel=L"z", ylabel=L"Normalized \chi", size=(800,600); plot_style...)
    plot!(z, D_growth_array ./ maximum(D_growth_array), label=L"D(z)", xlabel=L"z", ylabel=L"Normalized D(z)", title=L"Normalized D(z)", legend=:bottomright; plot_style...)
    plot!(z, n_chi_norm ./ maximum(n_chi_norm), label=L"n(z)", xlabel=L"z", ylabel=L"Normalized n(z)", title=L"Normalized n(z)", size=(800,600); plot_style...)
    plot!(z, Hubble_array ./ maximum(Hubble_array), label=L"H(z)", xlabel=L"z", ylabel=L"Normalized H(z)", title=L"Normalized H(z)", size=(800,600); plot_style...)
    plot!(z, W_x ./ maximum(W_x), label=L"\frac{n(z) b(z) D(z) \chi(z)^2 H(z)}{c}", xlabel=L"z", ylabel="Normalized Kernel", title="Normalized Kernel", legend=:bottomright; plot_style...)
    savefig(joinpath(output_dir, "plots/W/Wg_norm.png"))
    bb.append_to_log(output_dir, "Normalized kernel plot saved.")

    bb.append_to_log(output_dir, "Done kernel calculation.")
    bb.append_to_log(output_dir, "Size of the complete prefactor array is $(length(W_x)).")

    plot(z, W_x, label=L"\frac{n(z) b(z) D(z) \chi(z)^2 H(z)}{c}", xlabel=L"z", ylabel="Unnormalized Kernel", title="Unnormalized Kernel", legend=:bottomright; plot_style...)
    savefig(joinpath(output_dir, "plots/W/Wg_nonnorm.png"))
    bb.append_to_log(output_dir, "Unnormalized kernel plot saved.")

    npzwrite(joinpath(output_dir, "quantities/W/W_x.npy"), W_x)
    bb.append_to_log(output_dir, "W_x.npy saved.")

    return W_x, bz_array, D_growth_array, Hubble_array, n_chi_norm
end

function compute_Wcheb(
    xmin::Number, 
    xmax::Number, 
    n_cheb::Int,
    z::AbstractVector, 
    x::AbstractVector, 
    bias::AbstractVector, 
    growth::AbstractVector, 
    Hubble_param::AbstractVector, 
    n_chi_norm::AbstractVector; 
    output_dir::AbstractString, 
    sorting::Bool
)
    append_to_log(output_dir, "=== compute_Wcheb ===")
    append_to_log(output_dir, "The compute_Wcheb function has been called.")
    append_to_log(output_dir, "Computing Chebyshev nodes (xmin=$xmin, xmax=$xmax, n_cheb=$n_cheb, sorting=$sorting)...")

    if sorting
        chi_cheb_nodes = reverse(bb.get_clencurt_grid(xmin, xmax, n_cheb))
    else 
        chi_cheb_nodes = bb.get_clencurt_grid(xmin, xmax, n_cheb)
    end
    append_to_log(output_dir, "Size of the Chebyshev grid is $(length(chi_cheb_nodes)).")

    append_to_log(output_dir, "Computing interpolated values on Chebyshev grid...")
    
    # Redshift interpolation
    z_interp = DataInterpolations.AkimaInterpolation(z, x, extrapolation=ExtrapolationType.Linear)
    z_cheb_nodes = z_interp.(chi_cheb_nodes) 
    append_to_log(output_dir, "Size of the interpolated redshift array is $(length(z_cheb_nodes)).")

    # Bias interpolation
    b_interp = DataInterpolations.AkimaInterpolation(bias, x, extrapolation=ExtrapolationType.Linear)
    b_cheb = b_interp.(chi_cheb_nodes)
    append_to_log(output_dir, "Size of the interpolated bias array is $(length(b_cheb)).")

    # Growth factor interpolation
    D_interp = DataInterpolations.AkimaInterpolation(growth, x, extrapolation=ExtrapolationType.Linear)
    D_cheb = D_interp.(chi_cheb_nodes)
    append_to_log(output_dir, "Size of the interpolated growth factor array is $(length(D_cheb)).")

    # Normalization interpolation
    nz_interps = DataInterpolations.AkimaInterpolation(n_chi_norm, x, extrapolation=ExtrapolationType.Linear)
    nz_cheb = nz_interps.(chi_cheb_nodes)
    append_to_log(output_dir, "Size of the interpolated normalization array is $(length(nz_cheb)).")

    # Hubble parameter interpolation
    H_interp = DataInterpolations.AkimaInterpolation(Hubble_param, x, extrapolation=ExtrapolationType.Linear)
    Hubble_cheb = H_interp.(chi_cheb_nodes)
    append_to_log(output_dir, "Size of the interpolated Hubble factor array is $(length(Hubble_cheb)).")

    # Final prefactor computation
    W_cheb = chi_cheb_nodes.^2 .* b_cheb .* D_cheb .* nz_cheb .* Hubble_cheb
    append_to_log(output_dir, "Size of the calculated prefactor array W(χ) on Chebyshev nodes is $(length(W_cheb)).")

    npzwrite(joinpath(output_dir, "quantities/W_cheb/W_cheb.npy"), W_cheb)
    append_to_log(output_dir, "W_cheb.npy saved.")

    return W_cheb, chi_cheb_nodes
end

function compute_c_cheb(W_vals::AbstractVector; output_dir::AbstractString, sorting::Bool = true)
    append_to_log(output_dir, "=== compute_c_cheb ===")
    append_to_log(output_dir, "Planning FFT transformation for Chebyshev coefficients (sorting=$sorting)...")

    plan = bb.plan_fft(W_vals; sorting = sorting)
    append_to_log(output_dir, "Size of the planned FFT wave vector is $(length(plan)).")

    append_to_log(output_dir, "Computing Chebyshev coefficients via fast transformation...")
    c_cheb = bb.fast_chebcoefs(W_vals, plan; sorting = sorting)
    append_to_log(output_dir, "Size of the calculated Chebyshev coefficients vector is $(length(c_cheb)).")

    npzwrite(joinpath(output_dir, "quantities/c_cheb/c_cheb.npy"), c_cheb)
    append_to_log(output_dir, "c_cheb.npy saved.")

    return c_cheb
end
