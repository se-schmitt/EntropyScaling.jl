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

# Model container for pure transport property models 
#TODO replace by proper `split_model`
struct PureModelContainer{P<:AbstractTransportProperty,M} <: AbstractTransportPropertyModel
    prop::P
    puremodels::M
end

"""
    PureModelContainer <: AbstractTransportPropertyModel

Model container for pure transport property models. 
Calculating mixture properties is not possible with this model.
The interface for calculating transport properties is `property_name(model, p, T, x=Z1)`, where `x`` contains either the mole fractions `x::Vector` (only one entry may be not zero) or an integer (`x ∈ 1:N_components`).

## Parameters 

- `models`: transport property models for the pure components

## Example
```julia
using EntropyScaling

pure_model_A = ConstantModel(Viscosity(), 0.001)
pure_model_B = ConstantModel(Viscosity(), 0.002)
eta_model = PureModelContainer(Viscosity(), [pure_model_A, pure_model_B])

viscosity(eta_model, nothing, nothing, [0,1])
viscosity(eta_model, nothing, nothing, 2)

viscosity(eta_model, nothing, nothing, [0.5,0.5]) # Error
```
"""
PureModelContainer

function PureModelContainer(prop, ::Type{𝕄}, components; kwargs...) where {𝕄<:AESM}
    # models = [init_model(MODEL, c, userlocations, verbose) for c in components]
    models = [𝕄(c) for c in components]
    return PureModelContainer(prop, models)
end
function PureModelContainer(prop, models::Vector{𝕄}, components; kwargs...) where {𝕄<:Union{Nothing,ATPM}}
    return PureModelContainer(prop, models)
end

for (fun, prop) in PROPERTY_FUNCTIONS
    @eval begin
        function $fun(model::PureModelContainer{<:$prop,TT}, p, T, z=Z1; kwargs...) where {TT}
            TYPE = Base.promote_eltype(TT,p,T)
            z isa Vector && sum(iszero.(z)) != 1 && error("`PureModelContainer` only works for pure components!")
            z isa Number && length(z) > length(model.puremodels) && error("Substance ID too large!")
            idx_model = (z isa Number) ? z : findfirst((!).(iszero.(z)))
            pure_i = model.puremodels[idx_model]
            return isnothing(pure_i) ? TYPE(NaN) : $fun(pure_i, p, T; kwargs...)
        end
    end
end

export ConstantModel, PureModelContainer