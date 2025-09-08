module EntropyScalingMakieExt

using EntropyScaling
using Makie

const ES = EntropyScaling
const MK = Makie

# Helper function for plotting with the Makie interface
function _makie_plot(fig::MK.Figure, ax::MK.Axis, model::ES.AESM, data::ES.ATPD; kwargs...)
    slims = get(kwargs, :slims, nothing)  
    exp_data, mod_data = ES.calc_plot_data(model, data; slims)    

    # Settings for colormap
    cprop = get(kwargs, :cprop, nothing)
    cb_label = cprop == :T ? "T / K" : 
               cprop == :p ? "p / Pa" : 
               cprop == :ϱ ? "ϱ / mol/m³" : "$cprop"

    # Plot data points
    marker = get(kwargs, :marker, :circle)
    markersize = get(kwargs, :markersize, 15)
    markercolor = get(kwargs, :markercolor, :dodgerblue)
    colormap = get(kwargs, :colormap, :viridis)
    if isnothing(cprop)
        scatter!(ax, exp_data...; markersize, color=markercolor, marker)
    else
        cdata = getfield(data, cprop)
        colorrange = get(kwargs, :colorrange, extrema(cdata))
        p = scatter!(ax, exp_data...; markersize, color=cdata, colormap, colorrange)
        Colorbar(fig[1, 2], p, label=cb_label)
    end

    # Plot model line
    linecolor = get(kwargs, :linecolor, :black)
    linewidth = get(kwargs, :linewidth, 2)
    linestyle = get(kwargs, :linestyle, :solid)
    label = get(kwargs, :label, nothing)
    lines!(ax, mod_data...; linewidth, color = linecolor, linestyle, label)

    # Axis styling
    if isempty(ax.xlabel.val)
        ax.xlabel = "sˢ"
    end
    if isempty(ax.ylabel.val)
        ax.ylabel = "$(ES.symbol(data.prop))ˢ"
    end
    !isnothing(slims) ? xlims!(ax, slims) : nothing
    
    return fig, ax
end

# Define plot methods for Makie
function MK.plot(model::ES.AESM, data::ES.ATPD; kwargs...)
    fig = Figure()
    yscale = get(kwargs, :yscale, (data.prop == ThermalConductivity()) ? identity : log10)
    ax = Axis(fig[1, 1]; yscale)
    fig, ax = _makie_plot(fig, ax, model, data; kwargs...)
    return fig
end

function MK.plot!(ax, model::ES.AESM, data::ES.ATPD; kwargs...)
    fig = ax.parent
    _makie_plot(fig, ax, model, data; kwargs...)
    return fig
end

# Define a standalone plot! method that works with the current axis
function MK.plot!(model::ES.AESM, data::ES.ATPD; kwargs...)
    ax = current_axis()
    if isnothing(ax)
        fig = current_figure()
        if isnothing(fig)
            fig = Figure()
        end
        ax = Axis(fig[1, 1])
    else
        fig = ax.parent
    end
    _makie_plot(fig, ax, model, data; kwargs...)
    return fig
end

end # module
