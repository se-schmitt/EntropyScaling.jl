export GCESModel, GCESParams

struct GCESParams{P,T} <: AbstractEntropyScalingParams
    A::Vector{T}
    B::Vector{T}
    C::Vector{T}
    D::Vector{T}
    m::Vector{T}
    CE_model::ChapmanEnskogModel{T,<:AbstractCollisionIntegralMethod}
    base::BaseParam{P}
end

struct GCESModel{E, P, G} <: AbstractEntropyScalingModel
    components::Vector{<:AbstractString}
    groups::G
    params::P
    eos::E
end

"""
    GCESModel

Group contribution model for the viscosity from Löterging-Lin and Gross (2015) [lotgering-lin_group_2015](@cite).
The model is based on the homosegmented GC-PCP-SAFT EOS.

# Parameters

- `A` - `B`: component-specific parameters calculated from the groups

`m` (segment parameter of molecular-based EOS) are additional internal parameters (not to be set at 
construction).

# Example 

```julia
using EntropyScaling, Clapeyron, GCIdentifier

component = get_groups_from_smiles("CCCO", gcPCPSAFTGroups)
model = GCESModel(HomogcPCPSAFT(component), [component])

η = viscosity(model, 0.1e6, 300.)
```
"""
GCESModel

@modelmethods GCESModel GCESParams

function GCESModel(eos, components::Vector{<:Tuple})
    N_comps = length(components)
    names = first.(components)
    groups = last.(components)
    groups_eos = deepcopy(groups)

    prop = Viscosity()

    # Catch special cases
    for groups_i in groups
        if groups_i == ["CH3" => 2]       # ethane
            groups_i[:] .= ["CH3_eth" => 2]
        end
    end
    
    Aₐ, Bₐ, Cₐ, Dₐ = load_gc_params(GCESModel, prop, groups)

    mₐ = eos.params.segment
    σₐ = eos.params.sigma
    σᵢ = get_sig(eos)
    εᵢ = kB .* get_eps(eos)
    Mw = get_Mw(eos)
    mᵢ = get_m(eos)
    
    Aᵢ = zeros(N_comps)
    Bᵢ = zeros(N_comps)
    Cᵢ = zeros(N_comps)
    Dᵢ = zeros(N_comps)

    for (i,(groups_i,groups_eos_i)) in enumerate(zip(groups,groups_eos))
        γ = 0.45
        Vᵢ = 0
        for ((group, count), group_eos) in zip(groups_i, first.(groups_eos_i))
            Aᵢ[i] += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 * Aₐ[group]
            Bᵢ[i] += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 * Bₐ[group]
            Cᵢ[i] += count * Cₐ[group]
            Dᵢ[i] += count * Dₐ[group]

            Vᵢ += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 
        end
        Aᵢ[i] += log(sqrt(inv(mᵢ[i])))
        Bᵢ[i] /= (Vᵢ^γ)
    end
    
    CE_model = ChapmanEnskogModel(names,σᵢ,εᵢ,Mw,collision_integral=KimMonroe())
    params = GCESParams(Aᵢ, Bᵢ, Cᵢ, Dᵢ, mᵢ, CE_model, BaseParam(prop, Mw))

    return GCESModel(names, groups, params, eos)
end


function scaling_model(param::GCESParams{<:AbstractViscosity,T}, sˢ, x=Z1) where T
    m_mix = _dot(param.m, x)
    mx = x.*param.m./m_mix
    A = _dot(param.A, x)
    B = _dot(param.B, mx)
    C = _dot(param.C, mx)
    D = _dot(param.D, mx)

    ηˢ = exp(A + B*sˢ + C*(sˢ^2) + D*(sˢ^3))
    return ηˢ
end


function scaling(param::GCESParams, eos, Yˢ, T, ϱ, s, z; inv=true)
    k = inv ? 1 : -1
    prop = transport_property(param)
    Y₀ = property_CE(prop, param.CE_model, T, z)
    return Yˢ*Y₀^k
end

function scaling_variable(param::GCESParams, s, x = Z1)
    m_mix = _dot(param.m, x) 
    return s / R / m_mix
end