module ConfigMod 

using OrderedCollections

Base.@kwdef struct RunConfig
    # comoving distance integration parameters
    N::Int = 2^15 + 1
    xmin::Float64 = 26.0
    xmax::Float64 = 7000.0
    # k-grid
    Nk::Int = 150
    Nkp::Int = 150
    Nkpp::Int = 150
    n_cheb::Int = 200
    # ell values
    ℓ::Vector{Int} = collect(range(2, 200, length=100))
    # to sort the k-grid or not
    sorting::Bool = false
    # maximum wavenumber
    kmax::Float64 = 200.0 / 13.0
    # output directory
    output_root::String = "out/runs"
end

function validate(config::RunConfig)
    config.N > 1 || error("N must be > 1")
    config.xmin > 0 || error("xmin must be > 0")
    config.xmax > config.xmin || error("xmax must be larger than xmin")
    config.Nk > 0 || error("Nk must be > 0")
    config.Nkp > 0 || error("Nkp must be > 0")
    config.Nkpp > 0 || error("Nkpp must be > 0")
    config.n_cheb > 0 || error("n_cheb must be > 0")
    !isempty(config.ℓ) || error("ℓ must not be empty")
    minimum(config.ℓ) >= 0 || error("ℓ must be non-negative")
    return config
end
# minimum wavenumber for the k-grid
kmin(config::RunConfig) = 2.5 / config.xmax
# comoving distance grid
x_grid(config::RunConfig) = LinRange(config.xmin, config.xmax, config.N)

function k_grids(config::RunConfig, Blast)
    k  = Blast.get_clencurt_grid(kmin(config), config.kmax, config.Nk)
    kp = Blast.get_clencurt_grid(kmin(config), config.kmax, config.Nkp)
    kpp = Blast.get_clencurt_grid(kmin(config), config.kmax, config.Nkpp)

    if config.sorting
        return sort(k), sort(kp), sort(kpp)
    else
        return k, kp, kpp
    end
end


function write_run_config(
    output_dir,
    config,
    x,
    z,
    zmin,
    zmax,
    k_grid,
    kp_grid,
    kpp_grid
)
    kmin_val = 2.5 / config.xmax
    params_run = OrderedDict(
        "number of points for the integration of the integrals of W, N" => config.N,
        "minimum comoving distance xmin" => config.xmin,
        "maximum comoving distance xmax" => config.xmax,
        "minimum redshift zmin" => zmin,
        "maximum redshift zmax" => zmax,
        "minimum wavenumber kmin" => kmin_val,
        "maximum wavenumber kmax" => config.kmax,
        "number of Chebyshev nodes n_cheb" => config.n_cheb,
        "number of ℓ values" => length(config.ℓ),
        "number of k points Nk" => config.Nk,
        "number of kp points Nkp" => config.Nkp,
        "number of kpp points Nkpp" => config.Nkpp,
        "sorting enabled" => config.sorting
    )

    open(joinpath(output_dir, "config.txt"), "w") do io
        println(io, "=== Run parameters ===")
        for (key, val) in params_run
            println(io, key, " => ", val)
        end

        println(io)
        println(io, "=== Additional info ===")
        println(io, "comoving distance array x has size: ", length(x))
        println(io, "redshift array z has size: ", length(z))
        println(io, "k_grid and kp_grid are ", k_grid == kp_grid ? "the same" : "different")
        println(io, "k_grid and kpp_grid are ", k_grid == kpp_grid ? "the same" : "different")
        println(io, "kp_grid and kpp_grid are ", kp_grid == kpp_grid ? "the same" : "different")
        println(io, "k_grid is ", config.sorting ? "sorted" : "not sorted")
    end
end

end
