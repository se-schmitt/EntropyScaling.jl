struct ConstantModel{P,T} <: AbstractTransportPropertyModel
    prop::P 
    value::T
end

"""
    ConstantModel <: AbstractTransportPropertyModel

Convenient model for a constant-value transport property model.

## Parameters 

- `prop`: transport property
- `value`: transport property value

## Example
```julia
using EntropyScaling

eta_model = ConstantModel(Viscosity(), 0.001)
viscosity(eta_model, nothing, nothing)

D_model = ConstantModel(SelfDiffusionCoefficient(), 1e-9)
self_diffusion_coefficient(D_model, nothing, nothing)
```
"""
ConstantModel

for (fun, prop) in PROPERTY_FUNCTIONS
    @eval begin
        function $fun(model::ConstantModel{<:$prop,TT}, p, T, z=Z1; kwargs...) where {TT}
            x = z ./ sum(z)
            length(model.value) != length(z) && error("Number of constant values does not agree with moles!")
            return dot(model.value, x)
        end
    end
end


export ConstantModel