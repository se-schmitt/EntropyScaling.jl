module EntropyScalingPlotsExt

using EntropyScaling
using Plots

const ES = EntropyScaling

# Helper function for plotting with the Plots.jl interface
function _plots_plot(p::Plots.Plot, model::ES.AESM, data; kwargs...)
    slims = get(kwargs, :slims, nothing)
    prop = isnothing(data) ? kwargs[:prop] : data.prop 
    exp_data, mod_data = ES.calc_plot_data(model, data; slims, prop)    

    # Settings for colormap
    cprop = get(kwargs, :cprop, nothing)
    cb_label = cprop == :T ? "T / K" : 
               cprop == :p ? "p / Pa" : 
               cprop == :ϱ ? "ϱ / mol/m³" : "$cprop"

    # Plot data points
    if !isnothing(exp_data[1])
        marker = get(kwargs, :marker, :circle)
        markersize = get(kwargs, :markersize, 5)
        
        if isnothing(cprop)
            markercolor = get(kwargs, :markercolor, :blue)
            scatter!(p, exp_data...; markersize, color=markercolor, marker, label="", msw=0)
        else
            cdata = getfield(data, cprop)
            colormap = get(kwargs, :colormap, :viridis)
            colorrange = get(kwargs, :colorrange, extrema(cdata))
            scatter!(p, exp_data...; markersize, zcolor=cdata, seriescolor=colormap, clims=colorrange, marker, label="", msw=0)
            plot!(p; colorbar_title=cb_label, colorbar_framestyle=:box)
        end
    end

    # Plot model line
    linecolor = get(kwargs, :linecolor, :black)
    linewidth = get(kwargs, :linewidth, 2)
    linestyle = get(kwargs, :linestyle, :solid)
    label = get(kwargs, :label, "")
    plot!(p, mod_data...; linewidth, color=linecolor, linestyle, label)

    # Axis styling
    if !haskey(p.attr, :xlabel) || isempty(p.attr[:xlabel])
        plot!(p, xlabel="sˢ")
    end
    if !haskey(p.attr, :ylabel) || isempty(p.attr[:ylabel])
        plot!(p, ylabel="$(ES.symbol(prop))ˢ")
    end
    if !isnothing(slims)
        plot!(p, xlims=slims)
    end
    
    return p
end

# Define plot methods for Plots.jl
function Plots.plot(model::ES.AESM, data; kwargs...)
    p = plot()
    if !isnothing(data)
        yscale = (data.prop == ES.ThermalConductivity()) ? :identity : :log10
        plot!(p, yscale=yscale)
    end
    _plots_plot(p, model, data; kwargs...)
    return p
end

function Plots.plot!(p::Plots.Plot, model::ES.AESM, data; kwargs...)
    if !isnothing(data)
        yscale = (data.prop == ES.ThermalConductivity()) ? :identity : :log10
        plot!(p, yscale=yscale)
    end
    _plots_plot(p, model, data; kwargs...)
    return p
end

# Define a standalone plot! method that works with the current plot
function Plots.plot!(model::ES.AESM, data; kwargs...)
    p = current()
    if !isnothing(data)
        yscale = (data.prop == ES.ThermalConductivity()) ? :identity : :log10
        plot!(p, yscale=yscale)
    end
    _plots_plot(p, model, data; kwargs...)
    return p
end

end # module
