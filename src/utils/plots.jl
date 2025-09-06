"""
    plot(model::AbstractEntropyScalingModel, dat::TransportPropertyData; slims=nothing, 
         cprop=nothing)

Plots the scaled transport property as function of the entropy scaling variable (e.g. the reduced entropy).
The entropy scaling model as well as the scaled transport property data are shown.

`slims` sets the range for the entropy scaling variable (used for the model calculation).
`cprop` controls the property used for coloring (e.g. `:T`, `:p`, or `ϱ`).
A general documentation for `plot` for controling the appearance can be found [here](https://docs.juliaplots.org/stable/).

This functionality is provided by the `Plots.jl` extension. To use it, make sure to:
1. Import Plots.jl: `using Plots`
2. The plotting functionality will be automatically loaded
"""
plot
plot() = nothing