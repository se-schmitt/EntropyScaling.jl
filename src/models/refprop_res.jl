export RefpropRESModel, RefpropRESParams

"""
    RefpropRESParams
"""
struct RefpropRESParams{P,T} <: AbstractEntropyScalingParams
    n::Matrix{T}
    ξ::Vector{T}
    CE_model::ChapmanEnskogModel
    base::BaseParam{P}
end

function RefpropRESParams(prop::AbstractTransportProperty, n::Matrix{T}, ξ::Vector{T}, σ::Vector{T}, 
                          ε::Vector{T}, Mw::Vector{T}) where T
    
    CE_model = ChapmanEnskogModel(repeat([""],length(Mw)),σ,ε,Mw,collision_integral=KimMonroe())
    base = BaseParam(prop, Mw)
    return RefpropRESParams(n,ξ,CE_model,base)
end

"""
    RefpropRESModel{T} <: AbstractEntropyScalingModel

Entropy scaling model based on Refprop EOS.

## References
"""
struct RefpropRESModel{E,FP} <: AbstractEntropyScalingModel
    components::Vector{String}
    params::FP
    eos::E
end

@modelmethods RefpropRESModel RefpropRESParams

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

function scaling(param::RefpropRESParams, eos, Y, T, ϱ, s, z=Z1; inv=true)
    k = !inv ? 1 : -1
    Y₀ = property_CE(transport_property(param), param.CE_model, T, z)
    if inv
        return plus_scaling(param.base, Y, T, ϱ, s, z; inv=true) + Y₀
    else
        return plus_scaling(param.base, Y - Y₀, ϱ, s, z; inv=false)
    end
end
