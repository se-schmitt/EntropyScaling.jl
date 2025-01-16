"""
    plot(model::AbstractEntropyScalingModel, dat::TransportPropertyData; slims=nothing, 
         cprop=nothing)

Plots the scaled transport property as function of the entropy scaling variable (e.g. the reduced entropy).
The entropy scaling model as well as the scaled transport property data are shown.

`slims` sets the range for the entropy scaling variable (used for the model calculation).
`cprop` controls the property used for coloring (e.g. `:T`, `:p`, or `ϱ`).
A general documentation for `plot` for controling the appearance can be found [here](https://docs.juliaplots.org/stable/).
"""
plot
plot() = nothing

@recipe function f(model::AbstractEntropyScalingModel, dat::TransportPropertyData; 
                   slims=nothing, cprop=nothing )
    # Preprocessing 
    param = model[dat.prop]

    ϱdat = deepcopy(dat.ϱ)
    what_ϱ_nan = isnan.(ϱdat)
    ϱdat[what_ϱ_nan] = [molar_density(model.eos, dat.p[i], dat.T[i]; phase=dat.phase[i]) 
                        for i in findall(what_ϱ_nan)]
    sdat = entropy_conf.(model.eos, ϱdat, dat.T)
    sˢdat = scaling_variable.(param, sdat)
    Yˢdat = scaling.(param, model.eos, dat.Y, dat.T, ϱdat, sdat)

    slims = isnothing(slims) ? extrema(sˢdat) : slims
    sˢx = [range(slims...,100);]
    Yˢx = scaling_model.(param, sˢx)

    # Settings
    cb_title =  cprop == :T ? "T / K" : 
                cprop == :p ? "p / Pa" : 
                cprop == :ϱ ? "ϱ / mol/m³" : "$cprop"

    # General 
    yscale --> :log10
    framestyle --> :box
    xlabel --> "sˢ"
    ylabel --> "$(symbol(dat.prop))ˢ"
    xlims --> slims
    if !isnothing(cprop)
        colorbar_framestyle --> :box
        colorbar_title --> cb_title
    end

    # Plot data points 
    @series begin
        seriestype := :scatter
        msw --> 0
        markersize --> 5
        if isnothing(cprop)
            markercolor --> :dodgerblue
        else
            zcolor --> getfield(dat,cprop)
        end
        label := nothing
        z_order --> 1

        sˢdat, Yˢdat
    end

    # Plot model 
    @series begin
        seriestype := :line 
        linecolor --> :black
        linewidth --> 1.5
        label --> nothing

        sˢx, Yˢx
    end
    
    return nothing
end