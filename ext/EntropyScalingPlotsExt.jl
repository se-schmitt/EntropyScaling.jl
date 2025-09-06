module EntropyScalingPlotsExt

using EntropyScaling
import Plots

const ES = EntropyScaling

function Plots.plot(model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; 
                    slims=nothing, cprop=nothing, kwargs...)
    # Set default attributes for the initial plot
    plot_kwargs = Dict{Symbol, Any}()
    
    # Use log10 scale by default unless overridden
    if !haskey(kwargs, :yscale)
        plot_kwargs[:yscale] = :log10
    end
    
    # Set other default values if not specified
    if !haskey(kwargs, :framestyle)
        plot_kwargs[:framestyle] = :box
    end
    
    if !haskey(kwargs, :xlabel)
        plot_kwargs[:xlabel] = "sˢ"
    end
    
    if !haskey(kwargs, :ylabel)
        plot_kwargs[:ylabel] = "$(ES.symbol(dat.prop))ˢ"
    end
    
    if !haskey(kwargs, :xlims) && !isnothing(slims)
        plot_kwargs[:xlims] = slims
    end
    
    # Create a new plot with merged kwargs (user kwargs override defaults)
    p = Plots.plot(; merge(plot_kwargs, Dict(kwargs))...)
    
    # Add the model to the plot
    Plots.plot!(p, model, dat; slims=slims, cprop=cprop, kwargs...)
    
    return p
end

function Plots.plot!(p::Plots.Plot, model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; 
                    slims=nothing, cprop=nothing, kwargs...)
    # Preprocessing 
    param = model[dat.prop]

    ϱdat = deepcopy(dat.ϱ)
    what_ϱ_nan = isnan.(ϱdat)
    ϱdat[what_ϱ_nan] = [ES.molar_density(model.eos, dat.p[i], dat.T[i]; phase=dat.phase[i]) 
                        for i in findall(what_ϱ_nan)]
    sdat = ES.entropy_conf.(model.eos, ϱdat, dat.T)
    sˢdat = ES.scaling_variable.(param, sdat)
    Yˢdat = ES.scaling.(param, model.eos, dat.Y, dat.T, ϱdat, sdat)

    slims = isnothing(slims) ? extrema(sˢdat) : slims
    sˢx = [range(slims..., length=100);]
    Yˢx = ES.scaling_model.(param, sˢx)

    # Settings
    cb_title = cprop == :T ? "T / K" : 
               cprop == :p ? "p / Pa" : 
               cprop == :ϱ ? "ϱ / mol/m³" : "$cprop"
    
    # Add colorbar settings
    if !isnothing(cprop) && !haskey(kwargs, :colorbar_title)
        Plots.plot!(p; colorbar_title=cb_title, colorbar_framestyle=:box)
    end

    # Plot data points
    if isnothing(cprop)
        Plots.scatter!(p, sˢdat, Yˢdat;
            msw=0,
            markersize=5,
            label=get(kwargs, :label, nothing),
            markercolor=get(kwargs, :markercolor, :dodgerblue),
            z_order=1
        )
    else
        Plots.scatter!(p, sˢdat, Yˢdat;
            msw=0,
            markersize=5,
            label=get(kwargs, :label, nothing),
            zcolor=getfield(dat, cprop),
            z_order=1
        )
    end

    # Plot model line with either default black or the color specified in kwargs
    linecolor = get(kwargs, :lc, get(kwargs, :linecolor, :black))
    Plots.plot!(p, sˢx, Yˢx;
        linecolor=linecolor,
        linewidth=get(kwargs, :linewidth, 1.5),
        label=get(kwargs, :label, nothing)
    )
    
    return p
end

# Standalone plot! method that adds to the current plot
function Plots.plot!(model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; 
                     slims=nothing, cprop=nothing, kwargs...)
    # Get the current plot
    p = Plots.current()
    
    # Add to the current plot
    Plots.plot!(p, model, dat; slims=slims, cprop=cprop, kwargs...)
    
    return p
end

end # module