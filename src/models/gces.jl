export GCESModel, GCESParams

struct GCESParams{P,T} <: AbstractEntropyScalingParams
    mgcparams::Vector # mixing group-contribution parameters
    mᵢ::Vector{T}
    CE_model::ChapmanEnskogModel{T,<:AbstractCollisionIntegralMethod}
    base::BaseParam{P}
end

function GCESParams(prop::AbstractTransportProperty, eos, component_groups, mgcparams::Vector{Vector{Any}}, σᵢ , ϵᵢ, mᵢ, Mw)

    CE_model = ChapmanEnskogModel(first.(collect(component_groups)),σᵢ,ϵᵢ,Mw,collision_integral=KimMonroe())
    base = BaseParam(prop, Mw)
    return GCESParams(mgcparams, mᵢ, CE_model, base)
end 

"""
    GCESModel <: AbstractEntropyScalingModel

Group-contribution entropy scaling model [lotgering-lin_pure_2018](@cite) based on the homosegmented PCP-SAFT EOS [sauer_comparison_2014](@cite).
"""
struct GCESModel{E, P, H} <: AbstractEntropyScalingModel    #TODO adapt structure from Clapeyron.jl (same for params)
    components::Vector{<:AbstractString}
    groups::H
    params::P
    eos::E
end

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
    
    A = zeros(N_comps)
    B = zeros(N_comps)
    C = zeros(N_comps)
    D = zeros(N_comps)

    for (i,(groups_i,groups_eos_i)) in enumerate(zip(groups,groups_eos))
        γ = 0.45
        V_i = 0
        for ((group, count), group_eos) in zip(groups_i, first.(groups_eos_i))
            A[i] += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 * Aₐ[group]
            B[i] += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 * Bₐ[group]
            C[i] += count * Cₐ[group]
            D[i] += count * Dₐ[group]

            V_i += count * mₐ[group_eos] * (σₐ[group_eos] .* 1e10)^3 
        end
        A[i] += log(sqrt(inv(mᵢ[i])))
        B[i] /= (V_i^γ)
    end
    
    mgcparams = [A, B, C, D]
    CE_model = ChapmanEnskogModel(names,σᵢ,εᵢ,Mw,collision_integral=KimMonroe())
    params = GCESParams(mgcparams, mᵢ, CE_model, BaseParam(prop, Mw))

    return GCESModel(names, groups, params, eos)
end


function scaling_model(param::GCESParams{<:AbstractViscosity}, sˢ, x=Z1)
  
    A = param.mgcparams[1]
    B = param.mgcparams[2]
    C = param.mgcparams[3]
    D = param.mgcparams[4]

    ###TODO: hier fehlt formel um aus A_i, usw. in Vektoren richtige A_i, usw. zu machen, welche 
    ### dann in folgende formel eingesetzt werden.

    # Vorübergehend:
    A_i = A[1]
    B_i = B[1]
    C_i = C[1]
    D_i = D[1]
    ################

    ηˢ = exp(A_i + B_i*sˢ + C_i*(sˢ^2) + D_i*(sˢ^3))
    return ηˢ
end


function scaling(param::GCESParams, eos, ηˢ, T, ϱ, s, z; inv=true)

    # mᵢ = eos.pcpmodel.params.segment.values[1]
    
    ηₒ = viscosity(param.CE_model, T)
    # ηₒ = viscosity(param.CE_model, T) * sqrt(1/mᵢ)
    return ηˢ*ηₒ
end

function scaling_variable(param::GCESParams, s, z = Z1)
    sˢ = s/(NA*kB*param.mᵢ[1]) # m_gc entspricht m_gc,i aus dem Paper und sollte ein skalar sein, kein Vektor
    return sˢ
end