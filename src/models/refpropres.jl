export RefpropRES

struct CritTCNParam{T} <: AbstractParam
    φ0::CL.SingleParam{T}
    Γ::CL.SingleParam{T}
    qD::CL.SingleParam{T}
    Tref::CL.SingleParam{T}
end

abstract type AbstractRefpropRESParam{P,T} <: AbstractEntropyScalingParams{P} end

struct RefpropRESParam{P,T} <: AbstractRefpropRESParam{P,T}
    n1::CL.SingleParam{T}
    n2::CL.SingleParam{T}
    n3::CL.SingleParam{T}
    n4::CL.SingleParam{T}
    ξ::CL.SingleParam{T}
    crit::CritTCNParam{T}
    ce::ChapmanEnskog
    prop::P
end

struct RefpropRES{E,P} <: RefpropRESModel
    components::Vector{String}
    params::P
    eos::E
    sources::Vector{String}
end

# @modelmethods RefpropRES RefpropRESParam

"""
    RefpropRES{E,P} <: RefpropRESModel

    RefpropRES(components, eos=nothing; userlocations=String[], ce_userlocations=String[], verbose=false)

Entropy scaling model based on Refprop EOS for the viscosity and thermal conductivity [yang_linking_2022,yang_entropy_2021](@cite).

# Parameters

- `n1::SingleParam`: component-specific *or* global (group) parameters
- `n2::SingleParam`: component-specific *or* global (group) parameters
- `n3::SingleParam`: component-specific *or* global (group) parameters
- `n4::SingleParam`: component-specific *or* global (group) parameters
- `ξ::SingleParam`: component-specific scaling parameter (`ξ = 1` for individual fits)
- `φ0::CL.SingleParam`: critical parameter (optional, only for thermal conductivity)
- `Γ::CL.SingleParam`: critical parameter (optional, only for thermal conductivity)
- `qD::CL.SingleParam`: critical parameter (optional, only for thermal conductivity)
- `Tref::CL.SingleParam`: critical parameter (optional, only for thermal conductivity)

!!! info
    The default CoolProp EOS is used here which does not necessarily match the choice of the
    original papers. This might lead to slight deviations from values in the original papers
    (especially for the thermal conductivity).

# Example

```julia
using EntropyScaling

model_pure = RefpropRES("R134a")
η_pure = viscosity(model_pure, 1e5, 300.; phase=:liquid)

model_mix = RefpropRES(["decane","butane"])
η_mix = viscosity(model_mix, 1e5, 300., [.5,.5])
```
"""
RefpropRES

db_prefix(::Type{RefpropRES}) = "RefpropRES"
const PARAMS_REFPROPRES = ["ξ","n1","n2","n3","n4","φ0","Γ","qD","Tref"]
const CE_YANG2021_PATH = joinpath(DB_PATH,"ChapmanEnskogYang2021.csv")

function RefpropRES(components, eos=nothing; userlocations=String[], ce_userlocations=String[], verbose=false)
    _components = CL.format_components(components)

    _eos = _build_multifluid(_components, eos)
    
    params = RefpropRESParam[]
    for prop in [Viscosity(), ThermalConductivity()]
        _userlocations = prop in keys(userlocations) ? userlocations[prop] : String[]
        _params = CL.getparams(_components, [get_db_path(RefpropRES, prop)]; userlocations=_userlocations, ignore_missing_singleparams=PARAMS_REFPROPRES)
        components_missing = [all(_v.ismissingvalues[i] for (_,_v) in _params) for i in eachindex(_components)]

        if any(components_missing)
            verbose && @info "No RefpropRES $(name(prop)) parameters found for components: $(join(_components[components_missing],','))."
        else
            n1 = _params["n1"]
            n2 = _params["n2"]
            n3 = _params["n3"]
            n4 = _params["n4"]
            ξ = _params["ξ"]

            _ce_userlocations = isempty(ce_userlocations) ? CE_YANG2021_PATH : ce_userlocations 
            ce = ChapmanEnskog(_components; userlocations=_ce_userlocations, collision_integral=KimMonroe())
            crit = CritTCNParam(_params["φ0"], _params["Γ"], _params["qD"], _params["Tref"])

            push!(params, RefpropRESParam(n1,n2,n3,n4,ξ,crit,ce,prop))
        end
    end

    isempty(params) && error("No parameters found for components: $(join(components, ',')).")

    ref = ["10.1007/s10765-022-03096-9","10.1021/acs.iecr.1c02154"]

    return RefpropRES(_components, params, _eos, ref)
end

# Internal helper: build a MultiFluid EOS from component names
_build_multifluid(components, eos::CL.EoSModel) = eos
function _build_multifluid(components, ::Nothing)
    mixing = CL.init_model(CL.AsymmetricMixing, components, String[], false)
    eos = CL.MultiFluid(components; mixing)
    return eos
end

function scaling_model(param::RefpropRESParam, s, x=Z1)
    return exp(powerseries_scaling_model(param, s, x, (1.0, 1.5, 2.0, 2.5))) - 1.0
end

function scaling_model(param::RefpropRESParam{<:ThermalConductivity}, s, x=Z1)
    return powerseries_scaling_model(param, s, x, (1.0, 1.5, 2.0, 2.5))
end

function powerseries_scaling_model(param, s, x, g)
    Y⁺ = zero(Base.promote_eltype(s, x, g))
    for j in eachindex(x)
        Y⁺ += x[j] * (
            param.n1[j]*(s/param.ξ[j])^g[1] + 
            param.n2[j]*(s/param.ξ[j])^g[2] + 
            param.n3[j]*(s/param.ξ[j])^g[3] + 
            param.n4[j]*(s/param.ξ[j])^g[4]
        )
    end
    return Y⁺
end

function scaling(param::AbstractRefpropRESParam{P,TT}, eos, Y, T, ϱ, s, z=Z1; inv=true, η=nothing) where {P,TT}
    k  = !inv ? 1 : -1
    tp = transport_property(param)

    if tp == ThermalConductivity()
        each_ind = eachindex(z)
        λ₀   = [thermal_conductivity(param.ce, T; i) for i in each_ind]
        η₀   = [viscosity(param.ce, T; i) for i in each_ind]
        pure_eos = split_model(eos)
        cₚ₀  = isobaric_heat_capacity.(pure_eos, 1e-10, T)
        λ_int = thermal_conductivity_internal.(η₀, cₚ₀, get_Mw(eos))
        Δλ_c  = thermal_conductivity_critical(param.crit, eos, ϱ, T, η, z)
        Y₀    = mix_CE(MasonSaxena(), param.ce, λ₀ .+ λ_int, z; YΦ=λ₀) + Δλ_c
    else
        Y₀ = property_CE(tp, param.ce, T, z)
    end

    base = BaseParam(P(),param.ce.Mw)

    if inv
        return plus_scaling(base, Y, T, ϱ, s, z; inv=true) + Y₀
    else
        return plus_scaling(base, Y - Y₀, ϱ, s, z; inv=false)
    end
end

function ϱT_thermal_conductivity(model::RefpropRES, ϱ, T, z::AbstractVector=Z1)
    param = model[ThermalConductivity()]
    s  = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)
    λˢ = scaling_model(param, sˢ, z)
    η  = ϱT_viscosity(model, ϱ, T, z)
    return scaling(param, model.eos, λˢ, T, ϱ, s, z; inv=true, η)
end

function thermal_conductivity_internal(η₀, cₚ, Mw)
    f_int = 1.32
    return f_int * η₀/Mw * (cₚ - 5/2*R)
end

function thermal_conductivity_critical(crit, eos, ϱ, T, η, z)
    RD   = 1.02
    ν_γ  = 0.63 / 1.239

    _missing = any(any(getproperty(crit,k).ismissingvalues) for k in (:φ0, :Γ, :qD, :Tref))
    _missing && return 0.0

    φ0, Γ, qD, Tref = (_dot(getproperty(crit,k), z) for k in (:φ0, :Γ, :qD, :Tref))

    ∂ϱ∂p      = ForwardDiff.derivative(xp -> molar_density(eos, xp, T,    z; ϱ0=ϱ), pressure(eos, ϱ, T,    z))
    ∂ϱ∂p_Tref = ForwardDiff.derivative(xp -> molar_density(eos, xp, Tref, z; ϱ0=ϱ), pressure(eos, ϱ, Tref, z))
    Δ∂ϱ∂p = ∂ϱ∂p - Tref/T * ∂ϱ∂p_Tref
    if Δ∂ϱ∂p < 0
        return 0.0
    else
        cₚ     = isobaric_heat_capacity(eos, ϱ, T, z)
        cᵥ     = isochoric_heat_capacity(eos, ϱ, T, z)
        _, pc, ϱc = crit_mix(eos, z)

        cᵥ_cₚ  = cᵥ / cₚ
        ϱr     = ϱ / ϱc
        φ      = φ0 * (pc/ϱc * ϱr/Γ * Δ∂ϱ∂p)^(ν_γ)
        qDφ    = qD * φ
        ΔΩ     = 2/π * ((1-cᵥ_cₚ)*atan(qDφ) + cᵥ_cₚ*qDφ - 1 + exp(-qDφ/(1 + qDφ/3*(qDφ/ϱr)^2)))

        return ϱ*cₚ*RD*kB*T / (6π*η*φ) * ΔΩ
    end
end
