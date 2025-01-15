@recipe function f(model::AbstractEntropyScalingModel, dat::TransportPropertyData; 
                   slims=nothing, colprop=nothing )
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
    cb_title =  colprop == :T ? "T / K" : 
                colprop == :p ? "p / Pa" : 
                colprop == :ϱ ? "ϱ / mol/m³" : "$colprop"

    # General 
    yscale --> :log10
    framestyle --> :box
    xlabel --> "sˢ"
    ylabel --> "$(symbol(dat.prop))ˢ"
    xlims --> slims
    if !isnothing(colprop)
        colorbar_framestyle --> :box
        colorbar_title --> cb_title
    end

    # Plot data points 
    @series begin
        seriestype := :scatter
        msw --> 0
        markersize --> 5
        if isnothing(colprop)
            markercolor --> :dodgerblue
        else
            zcolor --> getfield(dat,colprop)
        end
        label --> nothing

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