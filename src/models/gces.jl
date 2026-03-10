using CSV, DataFrames

export GCESModel, GCESParams

struct GCESParams{P,T} <: AbstractEntropyScalingParams
    mgcparams::Vector # mixing group-contribution parameters
    CE_model::ChapmanEnskogModel{T,<:AbstractCollisionIntegralMethod}
    base::BaseParam{P}
end

function GCESParams(prop::AbstractTransportProperty, eos, n::Matrix{T}, m::Vector{T},
        σ::Vector{T}, ξ::Vector{T},ε::Vector{T}()) where T
    

    # Einlesen der VCParams
    df = CSV.read("C:/Users/leons/coding/Miniprojekt/MolekulareThermoMiniProjekt/database/viscosity_group_params.csv", DataFrame, header=false)

    A = []
    B = []
    C = []
    
    e = 1
    for ni in n
        i=1
        A_i = []
        B_i = []
        for nα in n[e]
            push!(A_i,nα* df[i, 2]) #dummy
            push!(B_i,nα*3* df[i, 3]) #dummy
            i+=1
            println(i)
        end
        push!(A, sum(A_i))
        push!(B, sum(B_i))
        e+=1
        
    end
    println(A,B)



    mgcparams = [[A, B], n, m , σ]
    Mw = convert(typeof(ξ),get_Mw(eos))
    CE_model = ChapmanEnskogModel(repeat([""],length(Mw)),σ,ε,Mw,collision_integral=KimMonroe())
    base = BaseParam(prop, Mw)
    return GCESParams(mgcparams, CE_model, base)
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

function GCESModel(eos, components::Vector{<:AbstractString}, groups)
    params =  GCESParams(prop, eos, n, m, σ, ξ, ε)
    return GCESModel(components, groups, params, eos)
end

# Define scaling methods
function scaling_model(param::GCESParams{<:AbstractViscosity}, s::TS, x=Z1)
    
end

function scaling(param::GCESParams, eos, Y, T, ϱ, s, z=Z1; inv=false)
    k = !inv ? 1 : -1

    # ...
end

function scaling_variable(param::GCESParams, s, z = Z1)
      # ...
end