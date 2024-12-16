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

get_m(m::RefpropRESParams) = FillArrays.Fill(1.0,length(m.ε))

"""
    RefpropRESModel{T} <: AbstractEntropyScalingModel

A generic entropy scaling model.
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
    return LogExpFunctions.logexpm1(lnY⁺p1) # Y⁺
end

function scaling(param::RefpropRESParams, eos, Y, T, ϱ, s, z=Z1; inv=true)
    k = !inv ? 1 : -1
    Y₀ = property_CE(param, T, z)
    if inv
        return plus_scaling(param.base, Y, T, ϱ, s, z; inv=true) + Y₀
    else
        return plus_scaling(param.base, Y - Y₀, ϱ, s, z; inv=false)
    end
end

viscosity_refprop(T,Mw,σ,ε) = viscosity_CE(T, Mw, σ, ε, Ω_22_neufeld(T*kB/ε))

function property_CE(param::RefpropRESParams{Viscosity}, T, z = Z1; mixing = nothing)
    Mw = param.base.Mw
    σ,ε = param.σ,param.ε
    prop = transport_property(param)
    if length(Mw) == 1
        Y₀ = viscosity_refprop(T,Mw[1],σ[1],ε[1])*one(eltype(z))
    else
        Y₀_all = viscosity_refprop.(T,Mw,σ,ε)
        Y₀ = mix_CE(mixing, param.base, Y₀⁺_all)
    end
    return Y₀
end
