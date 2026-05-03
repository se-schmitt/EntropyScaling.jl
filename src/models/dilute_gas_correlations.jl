export PolynomialDiluteGas

"""
    PolynomialDiluteGas <: AbstractDiluteGasModel

    PolynomialDiluteGas(components; userlocations=String[], prop=Viscosity())

Polynomial model for dilute gas transport properties (only temperature-dependent).

Currently, only parameters for the viscosity from [martinek_entropy_2025](@citet) are available in the database.

# Fields
- `m::NTuple{5,SingleParam{T}}`: polynomial coefficients mᵢ
- `Mw::SingleParam{T}`: molar mass (`[Mw] = g mol⁻¹`)

# Constructor

- `PolynomialDiluteGas(components; userlocations=String[], prop=Viscosity())`: database constructor

# Example

```julia
using EntropyScaling 

model = PolynomialDiluteGas(["butane","methanol"])
η_mix = viscosity(model, NaN, 300., [.5,.5])
```
"""
struct PolynomialDiluteGas{P,T} <: AbstractDiluteGasModel
    components::Vector{String}
    m::NTuple{5,CL.SingleParam{T}}
    Mw::CL.SingleParam{T}
    sources::Vector{String}
    prop::P
end

PolynomialDiluteGas

const POLYNOMIAL_DILUTE_GAS_COEFFICIENTS = ("m0", "m1", "m2", "m3", "m4")

db_model_path(::Type{PolynomialDiluteGas}) = joinpath("Others", "PolynomialDiluteGas_[PROP].csv")

# Constructor
function PolynomialDiluteGas(components; userlocations=String[], prop=Viscosity())
    _components = CL.format_components(components)

    params = CL.getparams(_components, [get_db_path(PolynomialDiluteGas, prop, nothing)]; userlocations)

    Mw = params["Mw"]
    m = ntuple(i -> params[POLYNOMIAL_DILUTE_GAS_COEFFICIENTS[i]], length(POLYNOMIAL_DILUTE_GAS_COEFFICIENTS))
    sources = hasproperty(Mw, :sourcecsvs) ? unique(String.(Mw.sourcecsvs)) : String[]

    return PolynomialDiluteGas(_components, m, Mw, sources, prop)
end
PolynomialDiluteGas(comps::String; kwargs...) = PolynomialDiluteGas([comps]; kwargs...)

function _property_polynomial(model::PolynomialDiluteGas{P,TT}, T, z) where {P,TT}
    x = z ./ sum(z)
    Y0_pure = zeros(promote_type(TT, typeof(T), eltype(x)), length(x))
    for (i, mi) in enumerate(model.m)
        Y0_pure .+= mi.values .* T.^(i - 1)
    end
    return x, Y0_pure
end

function viscosity(model::PolynomialDiluteGas{TP,TT}, p, T::T2, z::T3=Z1) where 
    {TP<:AbstractViscosity,TT,T2,T3}

    x, η0_pure = _property_polynomial(model, T, z)
    return mix_CE(Wilke(), model, η0_pure, x)
end

function thermal_conductivity(model::PolynomialDiluteGas{TP,TT}, p, T::T2, z::T3=Z1) where 
    {TP<:AbstractThermalConductivity,TT,T2,T3}

    x, λ0_pure = _property_polynomial(model, T, z)
    return mix_CE(MasonSaxena(), model, λ0_pure, x)
end
    