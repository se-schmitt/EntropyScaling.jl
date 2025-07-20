export PolynomialDiluteGasModel

"""
    PolynomialDiluteGasModel <: AbstractDiluteGasModel

Polynomial model for dilute gas tranasport properties (only temperature dependent).

Currently, only parameters for the viscosity from [martinek_entropy_2025](@citet) are available in the database.

# Fields
- `m::Matrix{T}`: polynomial coefficients mᵢ
- `Mw::Vector{T}`: molar mass (`[Mw] = kg mol⁻¹`)

# Constructor

- `PolynomialDiluteGasModel(components; ref="", ref_id="", prop=Viscosity())`: database constructor

# Example

```julia
using EntropyScaling 

model = PolynomialDiluteGasModel(["butane","methanol"])
η_mix = viscosity(model, NaN, 300., [.5,.5])
```
"""
struct PolynomialDiluteGasModel{P,T} <: AbstractDiluteGasModel
    components::Vector{String}
    m::Matrix{T}
    Mw::Vector{T}
    ref::Vector{Reference}
    prop::P
end

# Constructor
function PolynomialDiluteGasModel(comps::Vector{String}; ref="", ref_id="", prop=Viscosity())
    out = load_params(PolynomialDiluteGasModel, prop, comps; ref, ref_id)
    ismissing(out) ? throw(MissingException("No polynomial dilute gas parameters found for system [$(join(comps,", "))]")) : nothing
    Mw, _m..., refs = out 
    m = permutedims(hcat(_m...))
    return PolynomialDiluteGasModel(comps, m, Mw, refs, prop)
end
PolynomialDiluteGasModel(comps::String; kwargs...) = PolynomialDiluteGasModel([comps]; kwargs...)

function viscosity(model::PolynomialDiluteGasModel{TP,TT}, p, T::T2, z::T3=Z1) where 
    {TP<:AbstractViscosity,TT,T2,T3}

    x = z./sum(z)
    η0_pure = zeros(promote_type(TT,T2,eltype(T3)), length(x))
    for (i,mi) in enumerate(eachrow(model.m))
        η0_pure .+= mi .* T.^(i-1)
    end
    return mix_CE(Wilke(), model, η0_pure, x)
end

function thermal_conductivity(model::PolynomialDiluteGasModel{TP,TT}, p, T::T2, z::T3=Z1) where 
    {TP<:AbstractThermalConductivity,TT,T2,T3}

    x = z./sum(z)
    λ0_pure = zeros(promote_type(TT,T2,eltype(T3)), length(x))
    for (i,mi) in enumerate(eachrow(model.m))
        λ0_pure .+= mi .* T.^(i-1)
    end
    return mix_CE(MasonSaxena(), model, λ0_pure, x)
end
    