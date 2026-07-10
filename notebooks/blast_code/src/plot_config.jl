using Plots

"""
    setup_plot_theme() -> NamedTuple
"""
function setup_plot_theme(paper = nothing)
    # dimensions and resolution
    size_plot    = (1200, 600)
    size_heatmap = (600, 600)
    size_Cl = (1500, 800)
    if paper != nothing
        dpi          = 1000
    else
        dpi          = 300
    end
    
    # margins
    leftmargin   = 5Plots.mm
    bottommargin = 5Plots.mm
    rightmargin  = 5Plots.mm
    topmargin    = 5Plots.mm
    
    # Font delle etichette
    titlefontsize  = 15
    yguidefontsize = 15 
    xguidefontsize = 15
    legendfontsize = 10
    
    # Colori e palette
    c      = :curl
    colors = palette(:viridis)
    
    # Dizionario di stile riutilizzabile al volo nei plot: plot(...; shared_style...)
    shared_style = Dict{Symbol, Any}(
        :dpi            => dpi,
        :leftmargin     => leftmargin,
        :bottommargin   => bottommargin,
        :rightmargin    => rightmargin,
        :topmargin      => topmargin,
        :titlefontsize  => titlefontsize,
        :legendfontsize => legendfontsize
    )
    
    return (
        size_plot=size_plot, size_heatmap=size_heatmap, size_Cl=size_Cl, dpi=dpi,
        leftmargin=leftmargin, bottommargin=bottommargin, 
        rightmargin=rightmargin, topmargin=topmargin,
        titlefontsize=titlefontsize, yguidefontsize=yguidefontsize, 
        xguidefontsize=xguidefontsize, legendfontsize=legendfontsize,
        c=c, colors=colors, shared_style=shared_style
    )
end