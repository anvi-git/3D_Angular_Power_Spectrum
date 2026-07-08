module PathsMod

using Dates

Base.@kwdef struct RunPaths
    timestamp::String
    output_dir::String
    plot_subdir::String
    sl_plots::String
    kernel_plots::String
    quantity_subdir::String
    chebcoefs::String
    sl::String
end

function make_run_paths(; root="out/runs")
    timestamp = Dates.format(now(), "yyyy_mm_dd_HHMMSS")
    output_dir = joinpath(root, "run_$timestamp")
    plot_subdir = joinpath(output_dir, "plots")
    sl_plots = joinpath(plot_subdir, "Sl_plots")
    kernel_plots = joinpath(plot_subdir, "kernels")
    quantity_subdir = joinpath(output_dir, "quantities")
    chebcoefs = joinpath(quantity_subdir, "chebcoefs")
    sl = joinpath(quantity_subdir, "Sl")

    for path in (output_dir, plot_subdir, sl_plots, kernel_plots, quantity_subdir, chebcoefs, sl)
        mkpath(path)
    end

    return RunPaths(
        timestamp=timestamp,
        output_dir=output_dir,
        plot_subdir=plot_subdir,
        sl_plots=sl_plots,
        kernel_plots=kernel_plots,
        quantity_subdir=quantity_subdir,
        chebcoefs=chebcoefs,
        sl=sl,
    )
end

end
