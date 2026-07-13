using NPZ
using DataInterpolations

"""
    setup_cosmology_grid(; cosmo, N, Nk, Nkp, Nkpp, n_cheb) -> NamedTuple
"""
function setup_cosmology_grid(;
    cosmo  = Blast.FlatΛCDM(),
    N      = 2^15 + 1,
    Nk     = 350,
    Nkp    = 350,
    Nkpp   = 350,
    n_cheb = 200
)
    # background quantities
    z_b      = npzread("blast_code/data/background/z.npy")
    x_b      = npzread("blast_code/data/background/chi.npy")
    n5k_bins = npzread("blast_code/data/dNdzs_fullwidth.npz")
    
    # Akima interpolations
    z_of_χ   = DataInterpolations.AkimaInterpolation(z_b, x_b)
    chi_of_z = DataInterpolations.AkimaInterpolation(x_b, z_b)
    
    # grids in x and z
    xmin = 26
    xmax = 7000
    x    = LinRange(xmin, xmax, N)
    z    = z_of_χ.(x)
    
    zmin = minimum(z)
    zmax = maximum(z)
    
    # harmonic space and boundaries of the wavenumbers
    ℓ    = LinRange(2, 200, 100)
    kmax = 200 / 13 # 
    kmin = 2.5 / xmax
    
    return (
        cosmo=cosmo, N=N, Nk=Nk, Nkp=Nkp, Nkpp=Nkpp, n_cheb=n_cheb,
        z_b=z_b, x_b=x_b, n5k_bins=n5k_bins,
        z_of_χ=z_of_χ, chi_of_z=chi_of_z,
        xmin=xmin, xmax=xmax, x=x, z=z, zmin=zmin, zmax=zmax,
        ℓ=ℓ, kmax=kmax, kmin=kmin
    )
end

