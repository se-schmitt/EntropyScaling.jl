using PhysicalConstants.CODATA2018, EntropyScaling, Clapeyron

k = BoltzmannConstant

export GCESModel, GCESParams

struct GCESParams{P,T} <: AbstractEntropyScalingParams
    mgcparams::Vector # mixing group-contribution parameters
    mᵢ::Clapeyron.SingleParam{T}
    CE_model::ChapmanEnskogModel{T,<:AbstractCollisionIntegralMethod}
    base::BaseParam{P}
end

function GCESParams(prop::AbstractTransportProperty, eos, component_groups, mgcparams::Vector{Vector{Any}}, σᵢ , ϵᵢ, mᵢ, Mw) where T

    CE_model = ChapmanEnskogModel(first.(collect(component_groups)),σᵢ,ϵᵢ,Mw,collision_integral=KimMonroe())
    base = BaseParam(prop, Mw)
    return GCESParams(mgcparams, mᵢ, CE_model, base)
end 


"""
    GCESModel <: AbstractEntropyScalingModel

Group-contribution entropy scaling model [lotgering-lin_pure_2018](@cite) based on the homosegmented PCP-SAFT EOS [sauer_comparison_2014](@cite).
"""
struct GCESModel{E, P, H} <: AbstractEntropyScalingModel
    components::Vector{<:AbstractString}
    groups::H
    params::P
    eos::E
end

@modelmethods GCESModel GCESParams

function GCESModel(components, component_groups, eos_model)
    eos = eos_model.model

    ηref = nothing
    
    _ηref = isnothing(ηref) ? "Loetgering-Lin et al. (2015)" : ηref

    prop = Viscosity()
    
    out = load_params(GCESModel, prop, component_groups, GC = true)

    if !ismissing(out)
        A_a, B_a, C_a, D_a= out
    else
        throw(MissingException("No parameters found for system [$(join(components,", "))]"))
    end

    # base = BaseParam(prop, get_Mw(eos), refs)

    kB = EntropyScaling.kB

    m = eos.params.segment
    σ = eos.params.sigma
    σᵢ = [eos.pcpmodel.params.sigma.values[i,i] for i in 1:Int(sqrt(length(eos.pcpmodel.params.sigma.values)))]
    ϵᵢ = kB .* [eos.pcpmodel.params.epsilon.values[i,i] for i in 1:Int(sqrt(length(eos.pcpmodel.params.epsilon.values)))]
    mᵢ = eos.pcpmodel.params.segment
    Mw = eos.pcpmodel.params.Mw.values * 1e-3
    
    A = []  # In diesen Listen stehen später die Viskositätsparameter der Substanz i
    B = []  # hier sollen also die A_i, B_i, usw. drin stehen!
    C = []
    D = []
    
    for i in eachindex(component_groups) # somponents [i][1] wäre also substance i
        groups = component_groups[i][2] # groups ist dann die Liste der Gruppen in substance i
        
        A_i = 0
        B_i = 0
        C_i = 0
        D_i = 0
        V_i = 0

        γ = 0.45
        
        for group in eachindex(groups) # group = α (z.B. "-CH3") und count = n_α jeweils fur substance i
           
            A_i += groups[group][2] * m[group] * σ[group]^3 * A_a[group]
            B_i += groups[group][2] * m[group] * σ[group]^3 * B_a[group]
            C_i += groups[group][2] * C_a[group]
            D_i += groups[group][2] * D_a[group]

            V_i += groups[group][2] * m[group] * σ[group]^3 
        end
        B_i = B_i / (V_i^γ)

        push!(A,A_i)
        push!(B,B_i)
        push!(C,C_i)
        push!(D,D_i)
    end
    #######################################################################################
    
    mgcparams = [A, B, C, D]

    params = GCESParams(prop, eos, component_groups, mgcparams, σᵢ , ϵᵢ, mᵢ, Mw)

    return GCESModel(components, component_groups, params, eos)
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

    ηˢ = exp(A_i + B_i*sˢ + C_i*sˢ^2 + D_i*sˢ^3)
    return ηˢ
end


function scaling(param::GCESParams, _eos, ηˢ, T, ϱ, s, z; inv=true)
    ηₒ = viscosity(param.CE_model, T) * sqrt(param.mᵢ[1])###TODO: Warum Druck kein Unterschied?
    return ηˢ*ηₒ
end

function scaling_variable(param::GCESParams, s, z = Z1)
    sˢ = s/(NA*kB*param.mᵢ[1]) # m_gc entspricht m_gc,i aus dem Paper und sollte ein skalar sein, kein Vektor
    return sˢ
end