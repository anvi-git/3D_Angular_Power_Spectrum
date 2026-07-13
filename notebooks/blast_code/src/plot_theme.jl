module PlotThemeMod

using Plots

Base.@kwdef struct PlotTheme
    size_plot::Tuple{Int,Int} = (1200, 600)
    size_heatmap::Tuple{Int,Int} = (600, 600)
    dpi::Int = 300
    leftmargin::typeof(1Plots.mm) = 10Plots.mm
    bottommargin::typeof(1Plots.mm) = 10Plots.mm
    rightmargin::typeof(1Plots.mm) = 10Plots.mm
    topmargin::typeof(1Plots.mm) = 10Plots.mm
    titlefontsize::Int = 15
    yguidefontsize::Int = 15
    xguidefontsize::Int = 15
    legendfontsize::Int = 10
    colormap::Symbol = :curl
    palette_colors::Plots.ColorPalette
end

function make_plot_theme(n_colors::Int; palette_name::Symbol=:batlow)
    return PlotTheme(palette_colors = palette(palette_name, n_colors))
end

end
