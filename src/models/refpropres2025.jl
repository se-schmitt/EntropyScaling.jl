export RefpropRES2025

struct RefpropRES2025Param{P,T} <: AbstractRefpropRESParam{P,T}
    n1::CL.SingleParam{T}
    n2::CL.SingleParam{T}
    n3::CL.SingleParam{T}
    n4::CL.SingleParam{T}
    ξ::CL.SingleParam{T}
    ce::ChapmanEnskog
    prop::P
end

struct RefpropRES2025{E,P} <: RefpropRESModel
    components::Vector{String}
    params::P
    eos::E
    sources::Vector{String}
end

"""
    RefpropRES2025{E,P} <: RefpropRESModel

    RefpropRES2025(components, eos=nothing; userlocations=String[], ce_userlocations=String[], verbose=false)

Entropy scaling model based on Refprop EOS for the viscosity from Martinek et al. (2025) [martinek_entropy_2025](@cite) (see also [`RefpropRES`](@ref)).

# Parameters

- `n1::SingleParam`: component-specific *or* global (group) parameters
- `n2::SingleParam`: component-specific *or* global (group) parameters
- `n3::SingleParam`: component-specific *or* global (group) parameters
- `n4::SingleParam`: component-specific *or* global (group) parameters
- `ξ::SingleParam`: component-specific scaling parameter (`ξ = 1` for individual fits)

!!! info
    The default CoolProp EOS is used here which does not necessarily match the choice of the
    original papers. This might lead to slight deviations from values in the original papers
    (especially for the thermal conductivity).

# Example

```julia
using EntropyScaling

model_pure = RefpropRES2025("R134a")
η_pure = viscosity(model_pure, 1e5, 300.; phase=:liquid)

model_mix = RefpropRES2025(["decane","butane"])
η_mix = viscosity(model_mix, 1e5, 300., [.5,.5])
```
"""
RefpropRES2025

db_model_path(::Type{RefpropRES2025}) = joinpath("RefpropRES", "RefpropRES2025_[PROP].csv")

function RefpropRES2025(components, eos=nothing; userlocations=String[], ce_userlocations=String[], verbose=false)
    _components = CL.format_components(components)

    _eos = _build_multifluid(_components, eos)
    
    params = RefpropRES2025Param[]
    for prop in [Viscosity()]
        _userlocations = prop in keys(userlocations) ? userlocations[prop] : String[]
        _params = CL.getparams(_components, [get_db_path(RefpropRES2025, prop, nothing)]; userlocations=_userlocations, ignore_missing_singleparams=PARAMS_REFPROPRES)
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

            push!(params, RefpropRES2025Param(n1,n2,n3,n4,ξ,ce,prop))
        end
    end

    isempty(params) && error("No parameters found for components: $(join(components, ',')).")

    ref = ["10.1021/acs.jced.4c00451"]

    return RefpropRES(_components, params, _eos, ref)
end

function scaling_model(param::RefpropRES2025Param, s, x=Z1)
    return exp(powerseries_scaling_model(param, s, x, (1.8, 2.4, 2.8, 0.0))) - 1.0
end
