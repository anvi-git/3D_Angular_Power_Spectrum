using Dates
using OrderedCollections

function setup_output_directories()
    timestamp = Dates.format(now(), "yyyy_mm_dd_HHMMSS")
    
    # output directory structure
    output_dir       = "out/runs/run_$timestamp"
    plot_subdir      = joinpath(output_dir, "plots")
    Sl_plots         = joinpath(plot_subdir, "Sl_plots")
    Kernel_plots     = joinpath(plot_subdir, "kernels")
    quantity_subdir  = joinpath(output_dir, "quantities")
    chebcoefs        = joinpath(quantity_subdir, "chebcoefs")
    Sl               = joinpath(quantity_subdir, "Sl")
    
    # create directories
    mkpath(chebcoefs)
    mkpath(Sl)
    mkpath(Sl_plots)
    mkpath(Kernel_plots)
    
    println("Folders in place: ", output_dir, ", ", plot_subdir, " and ", quantity_subdir)
    
    # returns the paths
    return (
        output_dir      = output_dir,
        plot_subdir     = plot_subdir,
        Sl_plots        = Sl_plots,
        Kernel_plots    = Kernel_plots,
        quantity_subdir = quantity_subdir,
        chebcoefs       = chebcoefs,
        Sl              = Sl
    )
end

"""
    save_run_config(output_dir, N, xmin, xmax, zmin, zmax, kmin, kmax, n_cheb, ℓ, Nk, Nkp, Nkpp, x, z, k_grid, kp_grid, kpp_grid, sorting)

Crea un file `config.txt` nella cartella di output specificata, salvando tutti i metadati
e i parametri strutturali del run corrente.
"""
function save_run_config(
    output_dir, N, xmin, xmax, zmin, zmax, kmin, kmax, n_cheb, ℓ, Nk, Nkp, Nkpp, 
    x, z, k_grid, kp_grid, kpp_grid, sorting
)
    # Costruzione del dizionario ordinato con i parametri del run
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

    # Assicura che la cartella esista (anche se creata in precedenza)
    mkpath(output_dir)

    # Scrittura del file di testo
    open(joinpath(output_dir, "output_info.txt"), "w") do dictio
        println(dictio, "=== Run parameters ===")
        for (key, val) in params_run
            println(dictio, key, " => ", val)
        end
        println(dictio)
        println(dictio, "=== Additional info ===")
        println(dictio, "comoving distance array x has size: ", length(x))
        println(dictio, "redshift array z has size: ", length(z))

        # Controlli di uguaglianza tra le griglie
        println(dictio, "k_grid and kp_grid are ", k_grid != kp_grid ? "different" : "the same")
        println(dictio, "k_grid and kpp_grid are ", k_grid != kpp_grid ? "different" : "the same")
        println(dictio, "kp_grid and kpp_grid are ", kp_grid != kpp_grid ? "different" : "the same")

        # Controllo sull'ordinamento
        println(dictio, "k_grid is ", sorting ? "sorted" : "not sorted")
    end

    println("Data brought to safety in folder: ", output_dir)
    return params_run
end