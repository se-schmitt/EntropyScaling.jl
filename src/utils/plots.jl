"""
    plot(model::AbstractEntropyScalingModel, dat::TransportPropertyData; slims=nothing, 
         cprop=nothing, kwargs...)
    plot!(model::AbstractEntropyScalingModel, dat::TransportPropertyData; slims=nothing, 
          cprop=nothing, kwargs...)
    plot!(p::Plot, model::AbstractEntropyScalingModel, dat::TransportPropertyData; slims=nothing, 
          cprop=nothing, kwargs...)

Plots the scaled transport property as function of the entropy scaling variable (e.g. the reduced entropy).
The entropy scaling model as well as the scaled transport property data are shown.

## Arguments
- `model`: An entropy scaling model containing transport property parameters
- `dat`: Transport property data to plot
- `p`: (For `plot!` only) An existing plot to add the model to

## Keyword Arguments
- `slims`: Sets the range for the entropy scaling variable (used for the model calculation)
- `cprop`: Controls the property used for coloring (e.g. `:T`, `:p`, or `ϱ`)
- `kwargs`: Additional arguments passed to the plotting backend

## Available Functions
- `plot(model, dat)`: Creates a new plot with the model and data
- `plot!(model, dat)`: Adds the model and data to the current plot
- `plot!(p, model, dat)`: Adds the model and data to a specific plot

## Plot Backends
This functionality is provided by plotting backend extensions. Currently supported:
- **Plots.jl**: Import with `using Plots`
- **GLMakie.jl**: Import with `using GLMakie`

## Example
```julia
using EntropyScaling, Plots, Clapeyron

# Create a model and data
model = FrameworkModel(PCSAFT("propane"), [data])

# Create a plot with temperature coloring
plot(model, data; cprop=:T)

# Add a second model to the same plot with custom styling
model2 = FrameworkModel(PCSAFT("propane"), [data]; opts=FitOptions(...))
plot!(model2, data; lc=:blue, label="Model 2")
```
"""
plot
plot() = nothing

function calc_plot_data(model::AESM, data::TransportPropertyData; slims)
    param = model[data.prop]

    ϱdat = deepcopy(data.ϱ)
    what_ϱ_nan = isnan.(ϱdat)
    ϱdat[what_ϱ_nan] = [molar_density(model.eos, data.p[i], data.T[i]; phase=data.phase[i]) 
                        for i in findall(what_ϱ_nan)]
    sdat = entropy_conf.(model.eos, ϱdat, data.T)
    sˢdata = scaling_variable.(param, sdat)
    Yˢdata = scaling.(param, model.eos, data.Y, data.T, ϱdat, sdat)

    slims = isnothing(slims) ? extrema(sˢdata) : slims
    sˢx = [range(slims..., length=100);]
    Yˢx = scaling_model.(param, sˢx)
    
    return (sˢdata, Yˢdata), (sˢx, Yˢx)
end 