using Plots

"""
    setup_plot_theme() -> NamedTuple
"""
function setup_plot_theme(; paper = false)
    # dimensions and resolution
    size_plot    = (1200, 600)
    size_heatmap = (600, 600)
    size_Cl = (1500, 800)
    if paper == true
        dpi          = 1000
    else
        dpi          = 300
    end
    
    # margins
    leftmargin   = 10Plots.mm
    bottommargin = 10Plots.mm
    rightmargin  = 10Plots.mm
    topmargin    = 10Plots.mm
    
    # Font delle etichette
    titlefontsize  = 20
    yguidefontsize = 20 
    xguidefontsize = 20
    legendfontsize = 15
    
    # Colori e palette
    c      = :curl
    #colors = palette(:bamako)
    
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
        c=c, shared_style=shared_style
    )
end