"""
    plot(model::AbstractEntropyScalingModel, dat::TransportPropertyData; kwargs...)
    plot!(model::AbstractEntropyScalingModel, dat::TransportPropertyData; kwargs...)
    plot!(axis, model::AbstractEntropyScalingModel, dat::TransportPropertyData; kwargs...)

Plots the scaled transport property as function of the entropy scaling variable (e.g. the reduced entropy).
The entropy scaling model as well as the scaled transport property data are shown.

## Arguments
- `model`: An entropy scaling model containing transport property parameters
- `dat`: Transport property data to plot (or `nothing`)
- `axis`: An existing axis to add the plot

## Keyword Arguments
- `slims`: Sets the range for the entropy scaling variable (used for the model calculation)
- `prop`: Sets the transport property (only required if `dat == nothing`)
- `cprop`: Controls the property used for coloring (e.g. `:T`, `:p`, or `ϱ`)
- `marker`: Shape of markers for data points (default: `:circle`)
- `markersize`: Size of markers for data points (default: `15`)
- `markercolor`: Color of markers when not using property-based coloring (default: `:blue`)
- `colormap`: Colormap to use for coloring points (default: `:viridis`)
- `colorrange`: Range of values for the colorbar (default: extrema of the property used for coloring)
- `linecolor`: Color of the model line (default: `:black`)
- `linewidth`: Width of the model line (default: `2`)
- `linestyle`: Style of the model line (default: `:solid`)
- `label`: Label for legend entries (default: `nothing`)

## Supported Plot Packages
- [Makie.jl](https://docs.makie.org/stable/): Import a specific backend, e.g. `using GLMakie`.
- [Plots.jl](https://docs.juliaplots.org/stable/): Import with `using Plots`.

## Example
```julia
using EntropyScaling, Clapeyron, CoolProp
using GLMakie # or Plots

# Create synthetic data
sub = "propane"
T, p = rand(200.:500.,200), 10.0.^(rand(200).*4 .+ 4)
η = [PropsSI("V","T",T[i],"P",p[i],sub) for i in 1:200]

# Create data and model
ηdat = ViscosityData(T,p,nothing,η)
model_A = FrameworkModel(PCSAFT(sub), [ηdat])
model_B = FrameworkModel(PCSAFT(sub), [ηdat]; opts=FitOptions(what_fit=Dict(Viscosity() => Bool[0,1,0,1,0])))

# Plot 
fig = plot(model_A, ηdat; cprop=:T, linewidth=4, linecolor=:blue, label="model A (4 parameters)")
plot!(model_B, nothing; slims=(0,3), prop=Viscosity(), linestyle=:dash, label="model B (2 parameters)")
axislegend(position=:lt, framevisible=false)
ax2 = Axis(fig[2,1])
plot!(ax2, model_A, ηdat; marker=:star5, markercolor=:red, label="model A again")
axislegend(position=:lt, framevisible=false)
ax2.yscale = identity
fig
```
"""
plot
plot() = nothing

function calc_plot_data(model::AESM, data; slims, prop)
    param = model[prop]

    if !isnothing(data)
        ϱdat = deepcopy(data.ϱ)
        what_ϱ_nan = isnan.(ϱdat)
        ϱdat[what_ϱ_nan] = [molar_density(model.eos, data.p[i], data.T[i]; phase=data.phase[i]) 
                            for i in findall(what_ϱ_nan)]
        sdat = entropy_conf.(model.eos, ϱdat, data.T)
        sˢdata = scaling_variable.(param, sdat)
        Yˢdata = scaling.(param, model.eos, data.Y, data.T, ϱdat, sdat)
    else
        (sˢdata, Yˢdata) = (nothing, nothing)
    end
    
    @assert !(isnothing(slims) && isnothing(sˢdata)) "Please provide `slims` if `data` is `nothing`"
    slims = isnothing(slims) ? extrema(sˢdata) : slims
    sˢx = [range(slims..., length=100);]
    Yˢx = scaling_model.(param, sˢx)
    
    return (sˢdata, Yˢdata), (sˢx, Yˢx)
end 