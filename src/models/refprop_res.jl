export RefpropRESModel, RefpropRESParams

struct RefpropRESParams{P,T} <: AbstractEntropyScalingParams
    n::Matrix{T}
    ξ::Vector{T}
    CE_model::ChapmanEnskogModel
    base::BaseParam{P}
end

function RefpropRESParams(prop::AbstractTransportProperty, eos, n::Matrix{T}, ξ::Vector{T}, 
                          σ::Vector{T}, ε::Vector{T}) where T
    
    Mw = convert(typeof(ξ),get_Mw(eos))
    CE_model = ChapmanEnskogModel(repeat([""],length(Mw)),σ,ε,Mw,collision_integral=KimMonroe())
    base = BaseParam(prop, Mw)
    return RefpropRESParams(n,ξ,CE_model,base)
end

"""
    RefpropRESModel{T} <: AbstractEntropyScalingModel

Entropy scaling model based on Refprop EOS [yang_linking_2022](@cite). 

A database provides ready-to-use models for the viscosity of several fluids.
The model can favourably be used in combination with [`Clapeyron.jl`](https://github.com/ClapeyronThermo/Clapeyron.jl) and [`Coolprop.jl`](https://github.com/CoolProp/CoolProp.jl) (see examples).

## Parameters 

- `n::Matrix{T}`: component-specific *or* global (group) parameters
- `ξ::Vector{T}`: component-specific scaling parameter in case global parameters are used (`ξ = 1` for individual fits)
- `σ::Vector{T}`: LJ size parameter for the Chapman-Enskog model 
- `ε::Vector{T}`: LJ energy parameter for the Chapman-Enskog model

## Constructors

- `RefpropRESModel(eos, params::Dict{P})`: Default constructor (see above).
- `RefpropRESModel(eos, components)`: Creates a ES model using the parameters provided in the database (recommended). 
    `RefpropRESModel(components)` creates the EOS model on-the-fly (only works if `Clapeyron.jl` and `Coolprop.jl` are loaded).

## Example

```julia 
using EntropyScaling, Clapeyron, CoolProp

model_pure = RefpropRESModel("R134a")
η_pure = viscosity(model_pure, 1e5, 300.; phase=:liquid)

model_mix = RefpropRESModel(["decane","butane"])
η_mix = viscosity(model_mix, 1e5, 300., [.5,.5])
```
"""
struct RefpropRESModel{E,P} <: AbstractEntropyScalingModel
    components::Vector{String}
    params::P
    eos::E
end

@modelmethods RefpropRESModel RefpropRESParams

function RefpropRESModel(eos, components::Vector{String})
    params = RefpropRESParams[]
    for prop in [Viscosity()]   #TODO add thermal conductivity
        out = load_params(RefpropRESModel, prop, components)
        if !ismissing(out)
            ξ, n1, n2, n3, n4, refs = out
            CE_model = ChapmanEnskogModel(components; Mw=get_Mw(eos), ref_id="10.1007/s10765-022-03096-9")
            base = BaseParam(prop, get_Mw(eos), refs)
            push!(params, RefpropRESParams(permutedims(hcat(n1,n2,n3,n4)),ξ,CE_model,base))
        end
    end
    isempty(params) ? throw(MissingException("No parameters found for system [$(join(components,", "))]")) : nothing
    return RefpropRESModel(components, Tuple(params), eos)
end

function scaling_model(param::RefpropRESParams{Viscosity}, s, x=z1)
    g = (1.0,1.5,2.0,2.5)
    return generic_powerseries_scaling_model(param, s, x, g)
end

function generic_powerseries_scaling_model(param, s, x, g)
    n, ξ = param.n, param.ξ
    lnY⁺p1 = zero(Base.promote_eltype(n,s,x,g))
    @assert length(g) <= size(n,1)
    @assert length(x) == size(n,2)
    for i in eachindex(g)
        ni = zero(Base.promote_eltype(n,x))
        for j in eachindex(x)
            ni += x[j]*n[i,j]/(ξ[j]^g[i])
        end
        lnY⁺p1 += ni*s^g[i]
    end
    return LogExpFunctions.logexpm1(lnY⁺p1) # Y⁺
end

function scaling(param::RefpropRESParams, eos, Y, T, ϱ, s, z=Z1; inv=true, η_pure=nothing)
    k = !inv ? 1 : -1
    tp = transport_property(param)

    if tp == ThermalConductivity
        each_ind = eachindex(z)
        λ₀ = [thermal_conductivity(param.CE_model, T; i) for i in each_ind]
        η₀ = [viscosity(param.CE_model, T; i) for i in each_ind]
        pure_eos = split_model(eos)
        cₚ₀ = isobaric_heat_capacity.(pure_eos, 1e-10, T)
        λ_int = thermal_conductivity_internal.(η₀, cₚ, get_Mw(eos))
        cₚ = isobaric_heat_capacity.(pure_eos, ϱ, T)
        cᵥ = ischoric_heat_capacity.(pure_eos, ϱ, T)
        λ_c = thermal_conductivity_critical.()
        Y₀ = _Y₀ + λ_int + λ_c
    else
        Y₀ = property_CE(tp, param.CE_model, T, z)    
    end

    if inv
        return plus_scaling(param.base, Y, T, ϱ, s, z; inv=true) + Y₀
    else
        return plus_scaling(param.base, Y - Y₀, ϱ, s, z; inv=false)
    end
end

function ϱT_thermal_conductivity(model::RefpropRESModel, ϱ, T, z=Z1)
    param = model[ThermalConductivity()]
    s = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)
    λˢ = scaling_model(param, sˢ, z)
    η_pure = [ϱT_viscosity(model, ϱ, T, OneElement(length(z),i)) for i in eachindex(z)]
    return scaling(param, model.eos, λˢ, T, ϱ, s, z; inv=true, η_pure)
end

function thermal_conductivity_internal(η₀, cₚ, M)
    f_int = 1.32e-3
    return f_int*η₀/M * (cₚ - 5/2*R)
end

function thermal_conductivity_critical(ϱ, T, η, cₚ, cᵥ, )
    R₀ = 1.03
    Ω_c = 2.0/π * ((cₚ-cᵥ)/cₚ*atan(qd*ξ) + cᵥ/cₚ*qd*ξ)
    return ϱ*cₚ*R₀*kB*T/(6π*η*ξ) * (Ω_c-Ω₀_c)
end