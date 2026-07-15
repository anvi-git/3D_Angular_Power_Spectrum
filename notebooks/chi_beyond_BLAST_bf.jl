ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0;
using Pkg
Pkg.activate("blast_code"; io=devnull)
Pkg.resolve(; io=devnull)
Pkg.instantiate(; io=devnull)

using Revise 

using Base.Threads, NPZ, DataInterpolations, Interpolations, FastChebInterp
using BenchmarkTools, FFTW, FastTransforms, Dates, TOML, Plots, Plots.Measures
using QuadGK, LaTeXStrings, Tullio, StaticArrays, LoopVectorization, LinearAlgebra
using Unitful, SpecialFunctions, DifferentialEquations, Cosmology, NumericalIntegration, Trapz
using CSV, DataFrames, JSON, OrderedCollections
using MCIntegration
;

include("blast_code/src/Blast.jl")
include("blast_code/src/blast_tutorials.jl")
include("blast_code/src/galaxy_galaxy.jl")
include("blast_code/src/shear_shear.jl")
include("blast_code/src/config.jl")
include("blast_code/src/paths.jl")
include("blast_code/src/plot_config.jl")
using .galaxy_galaxy
using .shear_shear
using .Blast
using .blast_tutorials

import PhysicalConstants.CODATA2018: c_0 # speed of light in m/s

paths = setup_output_directories()
grid_data = setup_cosmology_grid() 
grids = Blast.generate_k_grids(grid_data.kmin, grid_data.kmax, grid_data.Nk, grid_data.Nkp, grid_data.Nkpp; sorting=false) 
plot_theme = setup_plot_theme() 

# ─────────────────────────────────────────────────────────────────────────────
# 1. COSMOLOGICAL FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
cosmo  = Blast.FlatΛCDM()
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s
ℓ    = LinRange(2, 200, 100)
N = 2^15+1
Nk = 150
Nkp = 150
Nkpp = 150
n_cheb = 200
# Comoving distance range supported by your redshift kernel n(x).
# Tune these to match your survey's x_min and x_max in Mpc/h or Mpc. 
xmin = 26.0   # Mpc/h
xmax = 7000.0  # Mpc/h
x    = LinRange(xmin, xmax, N)
##
z_b      = npzread("blast_code/data/background/z.npy")
x_b      = npzread("blast_code/data/background/chi.npy")
n5k_bins = npzread("blast_code/data/dNdzs_fullwidth.npz")
z_of_x   = DataInterpolations.AkimaInterpolation(z_b, x_b)
chi_of_z = DataInterpolations.AkimaInterpolation(x_b, z_b)
z    = z_of_x.(x)
kmax = 200 / 13
kmin = 2.5 / xmax
k = LinRange(kmin, kmax, Nk)
;

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
power_spectrum(k_pk, x1, x2) = @. sqrt(10^InterpPmm(z_of_x(x1),log10(k_pk)) * 10^InterpPmm(z_of_x(x2),log10(k_pk)))
power_spectrum_nl(k_pk, x1, x2) = @. sqrt(10^InterpPmm_nl(z_of_x(x1),log10(k_pk)) * 10^InterpPmm_nl(z_of_x(x2),log10(k_pk)))

Pk_grid = power_spectrum.(k, xmin, xmax);

# Hubble parameter H(x) in units where c=1 (or keep c explicit below)
Hubble_array = Blast.compute_hubble_factor.(z, Ref(cosmo)) ./ C_LIGHT;
#plot(z, Hubble_array, label = "Hubble factor", xlabel = L"z", ylabel = L"H(z) [Mpc/h]/c", legend = :topleft; plot_theme.shared_style...)

b_array = @. 0.98 + 1.24 * z - 1.72 * z^2 + 1.28 * z^3    ;
#plot(x, b_array, label = "Galaxy bias", xlabel = L"\chi [\mathrm{Mpc}/h]", ylabel = L"b(\chi)", legend = :topleft; plot_theme.shared_style...)

z_0 = 0.9/sqrt(2)
A = 1.5/z_0
alpha = 2
beta = 1.5
z_of_x = DataInterpolations.AkimaInterpolation(z, x); # z(x)
dz_dchi = [DataInterpolations.derivative(z_of_x, chi) for chi in x]
nz = A .* (z / z_0).^alpha .* exp.(-(z / z_0).^beta) # redshift distribution of sources, normalized to 1
nz_interp = DataInterpolations.AkimaInterpolation(nz, z, extrapolation=ExtrapolationType.Linear)
integral_nz, _ = quadgk(x -> nz_interp(x), minimum(z), maximum(z), rtol=1e-10)
nz_norm = nz ./ integral_nz
n_chi_norm = nz_norm .* dz_dchi 
;
#plot(z, n_chi_norm, label = "Normalized redshift distribution", xlabel = L"\chi [Mpc/h]", ylabel = L"n(\chi)", legend = :topleft; plot_theme.shared_style...)

D_growth_array = galaxy_galaxy.compute_growth_factor(cosmo, z);
#plot(z, D_growth_array, label = "Growth factor", xlabel = L"z", ylabel = L"D(z)", legend = :topleft; plot_theme.shared_style...)

Wchi = @. x^2 * Hubble_array * b_array * n_chi_norm * D_growth_array ;

gal_prefact_W, _, _, _, _ = galaxy_galaxy.galaxy_prefactor(grid_data.x, grid_data.z, grid_data.cosmo; 
                                                     output_dir=paths.output_dir, plot_style = plot_theme.shared_style);

Wchi_itp = DataInterpolations.AkimaInterpolation(Wchi, x)
;

println(grids.k_grid[1])
println(grids.k_grid[end])
krevgrid = @views reverse(grids.k_grid, dims=(1))
println(krevgrid[1])
println(krevgrid[end])


plot(x, gal_prefact_W)
plot!(x, Wchi)
plot!(x, Wchi_itp.(x))

println("Pk_grid is a vector of size $(size(Pk_grid))")


function Wtilde_quadgk(ℓ::Number, k::AbstractVector{<:Real}, k1::Real;
                               xmin=xmin, xmax=xmax,
                               rtol=1e-8, atol=0.0)

    Nk = length(k)
    vals = Vector{Float64}(undef, Nk)
    errs = Vector{Float64}(undef, Nk)

    Threads.@threads for i in eachindex(k)
        ki = k[i]

        fx(x) = Wchi_itp(x) *
                SpecialFunctions.sphericalbesselj(ℓ, ki * x) *
                SpecialFunctions.sphericalbesselj(ℓ, k1 * x)

        val, err = quadgk(fx, xmin, xmax; rtol=rtol, atol=atol)

        vals[i] = k1 * val
        errs[i] = k1 * err
        println("Wtilde_quadgk: ℓ=$ℓ, k1=$k1, ki=$ki, val=$(vals[i]), err=$(errs[i])")
    end

    return vals, errs
end

function S_l_bruteforce_quadgk(ℓ::Number, k1::Real, k2::Real, k::AbstractVector,
                               kmin_int, kmax_int,
                               xmin, xmax;
                               rtol_k = 1e-6, rtol_x = 1e-8,
                               atol_k = 0.0, atol_x = 0.0,
                               Nlg = 2/π)
    W1, _ = Wtilde_quadgk(ℓ, k, k1; xmin=xmin, xmax=xmax, rtol=rtol_x, atol=atol_x)
    #W2, _ = Wtilde_quadgk(ℓ, k, k2; xmin=xmin, xmax=xmax, rtol=rtol_x, atol=atol_x)
    
    # ∫ dk k^2 P(k) W1(k) W1(k)
    integrand = k.^2 .* Pk_grid .* W1 .* W1
    val = Nlg * trapz(k, integrand)
    
    return val, 0.0
end


function compute_S(ℓ::Number, k1s::AbstractVector, k2s::AbstractVector, 
                   k::AbstractVector,
                   kmin_int, kmax_int,
                   xmin, xmax;
                   rtol_k = 1e-6, rtol_x = 1e-8,
                   atol_k = 0.0, atol_x = 0.0,
                   Nlg = 2/π)
    Nk1 = length(k1s)
    Nk2 = length(k2s)
    S = Matrix{Float64}(undef, Nk1, Nk2)
    
    for (i1, k1) in enumerate(k1s)
        for (i2, k2) in enumerate(k2s)
            S[i1, i2], _ = S_l_bruteforce_quadgk(ℓ, k1, k2, k, kmin_int, kmax_int, 
                                                  xmin, xmax;
                                                  rtol_k=rtol_k, rtol_x=rtol_x,
                                                  atol_k=atol_k, atol_x=atol_x,
                                                  Nlg=Nlg)
        end
    end
    
    return S
end

# ─────────────────────────────────────────────────────────────────────────────
# Single-point benchmark: same k1 and k2
# ─────────────────────────────────────────────────────────────────────────────

ℓ_test = ℓ[2]
k_test = [0.01574015850603505]
k1_test = k_test
k2_test = k_test
println("k1_test = $k1_test h/Mpc, \nk2_test = $k2_test h/Mpc, \nℓ_test = $ℓ_test")

@time Sval = compute_S(
    ℓ_test, 
    k1_test, 
    k2_test, 
    krevgrid,
    kmin, kmax,
    xmin, xmax;
    rtol_k = 1e-6, 
    rtol_x = 1e-8,
    atol_k = 0.0, 
    atol_x = 0.0,
    Nlg = 2/π
)


println("Sval is a vector of FLoat64 of length ", length(Sval))
# S[k1 = 0.01574015850603505, k2 = 0.01574015850603505, l = 4.0] =  2.2330231701963857e-14
println(Sval)
println(size(Sval))
println("Sval[1] = ", Sval[1])


