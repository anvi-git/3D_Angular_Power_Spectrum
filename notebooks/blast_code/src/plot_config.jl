module PlotConfig

using Plots
using Plots.Measures
using LaTeXStrings

export SIZE_PLOT, SIZE_HEATMAP,
       PLOT_BASE, PLOT_PAPER, PLOT_PRESENTATION,
       HEATMAP_BASE, HEATMAP_PAPER,
       TWINX_BOTTOM, TWINX_TOP,
       apply_plot_paper!, apply_heatmap_paper!

const SIZE_PLOT = (1200, 600)
const SIZE_HEATMAP = (600, 600)

const PLOT_BASE = (
    size = SIZE_PLOT,
    dpi = 300,
    margin = 10mm,
    titlefontsize = 20,
)

const PLOT_PAPER = merge(PLOT_BASE, (
    xlabel = L"\ell",
    ylabel = L"S_\ell",
    titleposition = :left,
))

const PLOT_PRESENTATION = merge(PLOT_BASE, (
    size = (1000, 600),
    legendfontsize = 10,
    legendposition = :outertopright,
))

const HEATMAP_BASE = (
    size = SIZE_HEATMAP,
    dpi = 200,
    titlefontsize = 20,
    xguidefontsize = 15,
    yguidefontsize = 15,
)

const HEATMAP_PAPER = merge(HEATMAP_BASE, (
    xlabel = L"k",
    ylabel = L"k_p",
    xscale = :log10,
    yscale = :log10,
    c = :curl,
))

const TWINX_BOTTOM = (
    xaxis = L"\ell",
    ylabel = L"S_\ell",
    xscale = :log10,
)

const TWINX_TOP = (
    xaxis = L"k \; (h/\mathrm{Mpc})",
    xscale = :log10,
)

apply_plot_paper!(p; kws...) = plot!(p; PLOT_PAPER..., kws...)
apply_heatmap_paper!(args...; kws...) = heatmap(args...; HEATMAP_PAPER..., kws...)

end
