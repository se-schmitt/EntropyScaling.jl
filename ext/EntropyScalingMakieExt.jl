module EntropyScalingMakieExt

using EntropyScaling
import Makie

const ES = EntropyScaling

# Helper function for plotting with the Makie interface
function _makie_plot(
    model::ES.AbstractEntropyScalingModel, 
    dat::ES.TransportPropertyData;                
    slims=nothing, cprop=nothing, axis=nothing, figure=nothing, kwargs...
)                
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

    # Settings for colormap
    cb_title = cprop == :T ? "T / K" : 
               cprop == :p ? "p / Pa" : 
               cprop == :ϱ ? "ϱ / mol/m³" : "$cprop"

    # Create figure and axis if not provided
    fig = isnothing(figure) ? Makie.Figure() : figure
    
    if isnothing(axis)
        ax = Makie.Axis(fig[1, 1]; 
            xlabel = "sˢ",
            ylabel = "$(ES.symbol(dat.prop))ˢ",
            yscale = Makie.log10,
            limits = (slims, nothing)
        )
    else
        ax = axis
    end

    # Plot data points
    if isnothing(cprop)
        Makie.scatter!(ax, sˢdat, Yˢdat;
            markersize = 5,
            color = :dodgerblue,
            kwargs...)
    else
        cdata = getfield(dat, cprop)
        p = Makie.scatter!(ax, sˢdat, Yˢdat;
            markersize = 5,
            color = cdata,
            colormap = :viridis,
            kwargs...)
            
        # Add colorbar
        Makie.Colorbar(fig[1, 2], p, label=cb_title)
    end

    # Plot model line
    linecolor = get(kwargs, :linecolor, :black)
    Makie.lines!(ax, sˢx, Yˢx;
        linewidth = 2,
        color = linecolor)
    
    return fig, ax
end

# Define plot methods for Makie
function Makie.plot(model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; kwargs...)
    fig, ax = _makie_plot(model, dat; kwargs...)
    return fig
end

function Makie.plot!(axis, model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; kwargs...)
    _makie_plot(model, dat; axis=axis, figure=axis.parent, kwargs...)
    return axis.parent
end

# Define a standalone plot! method that works with the current axis
function Makie.plot!(model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; kwargs...)
    ax = Makie.current_axis()
    if isnothing(ax)
        fig = Makie.current_figure()
        if isnothing(fig)
            fig = Makie.Figure()
        end
        ax = Makie.Axis(fig[1, 1])
    end
    _makie_plot(model, dat; axis=ax, figure=ax.parent, kwargs...)
    return ax.parent
end

end # module
