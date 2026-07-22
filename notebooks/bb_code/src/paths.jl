function setup_output_directories()
    timestamp = Dates.format(now(), "yyyy_mm_dd_HHMMSS")
    
    # output directory structure
    out_dir       = "out/runs/run_$timestamp"
    plot_subdir      = joinpath(out_dir, "plots")
    W_plots          = joinpath(plot_subdir, "W")
    k_grid_plots     = joinpath(plot_subdir, "k_grid")
    Sl_plots         = joinpath(plot_subdir, "Sl_plots")
    heatmaps         = joinpath(Sl_plots, "heatmaps")
    quantity_subdir  = joinpath(out_dir, "quantities")
    c_cheb           = joinpath(quantity_subdir, "c_cheb")
    Sl               = joinpath(quantity_subdir, "Sl")
    W_cheb           = joinpath(quantity_subdir, "W_cheb")
    W                = joinpath(quantity_subdir, "W")
    
    # create directories
    mkpath(c_cheb)
    mkpath(Sl)
    mkpath(Sl_plots)
    mkpath(heatmaps)
    mkpath(W_plots)
    mkpath(k_grid_plots)
    mkpath(W_cheb)
    mkpath(W)

    println("Folders in place: \n", out_dir, "\n", plot_subdir, "\n and \n", quantity_subdir)
    
    # returns the paths
    return (
        output_dir      = out_dir,
        plot_subdir     = plot_subdir,
        Sl_plots        = Sl_plots,
        W_plots         = W_plots,
        k_grid_plots    = k_grid_plots,
        quantity_subdir = quantity_subdir,
        c_cheb          = c_cheb,
        Sl              = Sl,
        W_cheb          = W_cheb,
        W               = W,
        heatmaps        = heatmaps
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

    println("Data brought to safety in folder: \n", output_dir)
    append_to_log(output_dir, "Run started at $(Dates.now())")
    return params_run
end

"""
    append_to_log(output_dir, message)
"""
function append_to_log(output_dir, message::String)
    log_file = joinpath(output_dir, "output_info.txt")
    
    open(log_file, "a") do io
        time_str = Dates.format(now(), "HH:MM:SS")
        println(io, "[$time_str] ", message)
    end
end
