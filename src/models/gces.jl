export GCES

abstract type GCESModel <: AbstractEntropyScalingModel end

struct GCESParams{P,T} <: AbstractEntropyScalingParam{P}
    A::CL.SingleParam{T}
    B::CL.SingleParam{T}
    C::CL.SingleParam{T}
    D::CL.SingleParam{T}
    m::CL.SingleParam{T}
    ce::ChapmanEnskog{T,<:AbstractCollisionIntegralMethod}
    prop::P
end

struct GCES{E,P} <: GCESModel
    components::Vector{<:AbstractString}
    groups::CL.GroupParam{Float64}
    params::P
    eos::E
    sources::Vector{String}
end

"""
    GCESModel

Group contribution model for the viscosity [lotgering-lin_group_2015](@cite) and the thermal conductivity [hopp_thermal_2019](@cite).
The model is based on the homosegmented GC-PCP-SAFT EOS.

# Parameters

- `A` - `B`: component-specific parameters calculated from the groups

`m` (segment parameter of molecular-based EOS) are additional internal parameters (not to be set at 
construction).

# Example 

```julia
using EntropyScaling, Clapeyron, GCIdentifier

component = get_groups_from_smiles("CCCO", gcPCPSAFTGroups)
model = GCES(component, HomogcPCPSAFT(component))

η = viscosity(model, 0.1e6, 300.)
λ = thermal_conductivity(model, 0.1e6, 300.)
```
"""
GCES

db_model_path(::Type{GCES}) = joinpath("GCES", "GCES_[PROP].csv")
const PARAMS_GCES = ["A", "B", "C", "D"]

function GCES(components, eos=nothing; userlocations=Dict(), verbose=false)
    _components = CL.format_components(components)
    γ = 0.45

    _eos = _build_homogc_pcsaft(_components, eos)
    groups = _eos.groups
    groups_eos = deepcopy(groups)
    mₐ = _eos.params.segment
    σₐ = _eos.params.sigma
    σᵢ = _eos.pcpmodel.params.sigma
    εᵢ = _eos.pcpmodel.params.epsilon
    mᵢ = _eos.pcpmodel.params.segment

    # Catch special cases
    idx_eth = findfirst(groups.groups .== Ref(["CH3"]))
    if !isnothing(idx_eth)
        _groups = setindex!(copy(groups.groups), ["CH3_eth"], idx_eth)
        _N_groups_old = length(groups.flattenedgroups)
        _flattenedgroups = vcat(copy(groups.flattenedgroups), "CH3_eth")
        _n_flattenedgroups = [i == idx_eth ? vcat(zeros(_N_groups_old), 2.0) : vcat(v, 0.0) for (i,v) in enumerate(groups.n_flattenedgroups)]
        _i_groups = findall.(.!iszero.(v) for v in _n_flattenedgroups)
        groups = CL.GroupParam(_components, _groups, :gcPCPSAFT_ES, groups.n_groups, groups.n_intergroups, _i_groups, _flattenedgroups, _n_flattenedgroups, groups.sourcecsvs)
    end

    params_dict = _get_empty_params_dict()
    for prop in [Viscosity(), ThermalConductivity()]
        _userlocations = prop in keys(userlocations) ? userlocations[prop] : String[]
        _params = CL.getparams(groups, [get_db_path(GCES, prop, nothing)]; userlocations=_userlocations, ignore_missing_singleparams=PARAMS_GCES)
        components_missing = [all(_v.ismissingvalues[i] for (_, _v) in _params) for i in eachindex(_components)]

        if any(components_missing)
            verbose && @info "No RefpropRES $(name(prop)) parameters found for components: $(join(_components[components_missing],','))."
        else
            Aₐ = _params["A"]
            Bₐ = _params["B"]
            Cₐ = _params["C"]
            Dₐ = _params["D"]

            Aᵢ = CL.SingleParam("A", _components, zeros(length(_components)))
            Bᵢ = CL.SingleParam("B", _components, zeros(length(_components)))
            Cᵢ = CL.SingleParam("C", _components, zeros(length(_components)))
            Dᵢ = CL.SingleParam("D", _components, zeros(length(_components)))

            for (i,(n_i,grps_i,grps_eos_i)) in enumerate(zip(groups.n_groups, groups.groups, groups_eos.groups,))
                Vᵢ = 0
                n_total = sum(n_i)
                for (count, group, group_eos) in zip(n_i, grps_i, grps_eos_i)
                    if prop isa AbstractViscosity
                        Aᵢ[i] += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 * Aₐ[group]
                        Bᵢ[i] += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 * Bₐ[group]
                    elseif prop isa AbstractThermalConductivity
                        Aᵢ[i] += count * Aₐ[group]
                        Bᵢ[i] += count * Bₐ[group]
                    end

                    Cᵢ[i] += count * Cₐ[group]
                    if prop isa AbstractViscosity
                        Dᵢ[i] += count * Dₐ[group]
                    elseif prop isa AbstractThermalConductivity
                        Dᵢ[i] += n_total * Dₐ[group]
                    end
                    Vᵢ += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3
                end

                if prop isa AbstractViscosity
                    Aᵢ[i] += log(sqrt(inv(mᵢ[i])))
                    Bᵢ[i] /= (Vᵢ^γ)
                end
            end

            ce = ChapmanEnskog(_eos; collision_integral=KimMonroe())

            params_dict[prop] = GCESParams(Aᵢ,Bᵢ,Cᵢ,Dᵢ,mᵢ,ce,prop)
        end
    end

    params = ParamVector(params_dict)
    ismissing(params) && error("No parameters found for components: $(join(components, ',')).")

    ref = ["10.1021/acs.iecr.5b01698"]

    return GCES(_components, groups, params, _eos, ref)
end

# Internal helper: build a MultiFluid EOS from component names
_build_homogc_pcsaft(components, eos::CL.EoSModel) = eos
function _build_homogc_pcsaft(components, ::Nothing)
    eos = CL.HomogcPCPSAFT(components; assoc_options=CL.AssocOptions(; combining=:cr1))
    return eos
end

function scaling_model(param::GCESParams{<:AbstractViscosity,T}, sˢ, x=Z1) where {T}
    m_mix = _dot(param.m, x)
    mx = x .* param.m ./ m_mix
    A = _dot(param.A, x)
    B = _dot(param.B, mx)
    C = _dot(param.C, mx)
    D = _dot(param.D, mx)

    ηˢ = exp(A + B*sˢ + C*(sˢ^2) + D*(sˢ^3))
    return ηˢ
end

function scaling_model(param::GCESParams{<:AbstractThermalConductivity,T}, sˢ, x=Z1) where {T}
    m_mix = _dot(param.m, x)
    mx = x .* param.m ./ m_mix
    A = _dot(param.A, x)
    B = _dot(param.B, mx)
    C = _dot(param.C, mx)
    D = _dot(param.D, mx)

    λˢ = exp(A + B*sˢ + C*(1 - exp(sˢ)) + D*sˢ^2)
    return λˢ
end

function scaling(param::GCESParams, eos, Yˢ, T, ϱ, s, z; inverse=true)
    k = inverse ? 1 : -1
    prop = transport_property(param)
    Y₀ = property_CE(prop, param.ce, T, z)
    return Yˢ * Y₀^k
end

function scaling(param::GCESParams{<:AbstractThermalConductivity,F}, eos, Yˢ, T, ϱ, s, z; inverse=true) where {F}
    k = inverse ? 1 : -1
    prop = transport_property(param)

    x = z ./ sum(z)
    m_mix = _dot(param.m, x)
    σ_mix = sum(x[i] * eos.pcpmodel.params.sigma[i, i] for i in eachindex(x))
    ε_mix = sum(x[i] * eos.pcpmodel.params.epsilon[i, i] for i in eachindex(x))

    λ_CE = property_CE(prop, param.ce, T, z) * sqrt(m_mix)

    Tˢ = T / (ε_mix * m_mix)

    c1 = -0.0167141
    c2 = 0.0470581
    λ_int = (m_mix^2 * σ_mix^3 * ε_mix) * (c1 * Tˢ + c2 * Tˢ^2) * 1e25

    sˢ = scaling_variable(param, s, z)
    φ = exp(2 * sˢ)

    λ_ref = λ_CE + φ * λ_int

    return Yˢ * λ_ref^k
end

function scaling_variable(param::GCESParams, s, x=Z1)
    m_mix = _dot(param.m, x)
    return s / R / m_mix
end
