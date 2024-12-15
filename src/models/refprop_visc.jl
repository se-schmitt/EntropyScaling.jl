export RefpropRESModel, RefpropRESParams

"""
    RefpropRESParams
"""
struct RefpropRESParams{P,T} <: AbstractEntropyScalingParams
    n::Matrix{T}
    σ::Vector{T}
    ε::Vector{T}
    base::BaseParam{P}
end

"""
    FrameworkModel{T} <: AbstractEntropyScalingModel

A generic entropy scaling model.
"""
struct RefpropRESModel{E} <: AbstractEntropyScalingModel
    components::Vector{String}
    params::RefpropRESParams
    eos::E
end

function scaling_model(param::RefpropRESParams, s, x=[1.]) where T
    g = (1.0,1.5,2.0,2.5)
    return generic_powerseries_scaling_model(param, s, x, g)
end

property_CE(prop::Viscosity, T, Mw, σ, ε, Ω_22) = viscosity_CE(T, Mw, σ, ε, Ω_22)

function generic_powerseries_scaling_model(param::RefpropRESParams, s, x, g)
    ##α = a1,a2,a3,a4,ξ
    #ln(η + 1) = a1*(s/ξ)^g1 + a2*(s/ξ)^g2 + a3*(s/ξ)^g3 + a4*(s/ξ)^g4
    n = param.n
    lnY⁺p1 = zero(Base.promote_eltype(n,s,x,g))
    @assert length(g) <= size(n,1)
    @assert length(x) == size(n,2)
    for i in 1:length(g)
        ni = zero(Base.promote_eltype(n,x))
        for j in 1:length(x)
            ni += x[j]*n[i,j]
        end
        lnY⁺p1 += ni*s^g[i]
    end
    return exp(lnY⁺p1) - 1.         # Y⁺
end

function scaling(param::RefpropRESParams, eos, Y, T, ϱ, s, z=[1.]; inv=true)
    k = !inv ? 1 : -1
    
    if length(z) == 1
        _1 = one(eltype(z))
        Y₀ = _1*property_CE(transport_property(param.base), T, param.base.Mw[1], param.σ[1], param.ε[1], Ω_22_newfeld(T*kB/param.ε[1]))   #TODO make generic (only valid for viscosity currently)
    else
        Y₀_all = property_CE.(transport_property(param.base), T, param.base.Mw, param.σ, param.ε)
        Y₀ = mix_CE(param.base, η₀_all, z)
    end

    if inv
        return plus_scaling(param.base, Y, T, ϱ, s, z; inv=true) + Y₀
    else
        return plus_scaling(param.base, Y - Y₀, ϱ, s, z; inv=false)
    end
end

function ϱT_viscosity(model::RefpropRESModel, ϱ, T, z=[1.])
    s = entropy_conf(model.eos, ϱ, T, z)
    @show s
    η⁺p1 = scaling_model(model.params, -s/R, z)
    return scaling(model.params, model.eos, η⁺p1, T, ϱ, s, z)
end