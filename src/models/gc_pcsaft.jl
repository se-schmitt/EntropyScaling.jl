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

"""
    GCESModel <: AbstractEntropyScalingModel

Group-contribution entropy scaling model [lotgering-lin_pure_2018](@cite) based on the homosegmented PCP-SAFT EOS [sauer_comparison_2014](@cite).


"""
struct GCESModel{E,P} <: AbstractEntropyScalingModel
    components::Vector{<:AbstractString}
    params::P
    eos::E
end

@modelmethods GCESModel GCESParams

function GCESModel(eos)
    get_groups(eos)
    # ...
end

# Define scaling methods
function scaling_model(param::GCESParams{<:AbstractViscosity}, s::TS, x=Z1)
    sˢ = scaling_variable(param,s,x)
    m_ = _dot(param.m,z)
    P = zeros(TS,4)
    for i in eachindex(z)
        xᵢmᵢ_m = param.m[i]*z[i]/m_ 
        P[1] += z[i]*param.A[i]
        P[2] += xᵢmᵢ_m*param.B[i]
        P[3] += xᵢmᵢ_m*param.C[i]
        P[4] += xᵢmᵢ_m*param.D[i]
    end

    return exp(eval_poly(sˢ, P))
end

function scaling(param::GCESParams, eos, Y, T, ϱ, s, z=Z1; inv=false)
    k = !inv ? 1 : -1

    # Transport property scaling
    prop = transport_property(param)
    Y₀ = property_CE(prop, param.CE_model, T, z) / sqrt(get_m(param))
    
    # Entropy
    Yˢ = Y * Y₀^k

    return Yˢ
end

function scaling_variable(param::GCESParams, s, z = Z1)
    return -s / R / _dot(get_m(param),z)
end