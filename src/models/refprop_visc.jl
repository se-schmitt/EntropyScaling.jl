struct ViscosityREFPROP <: AbstractViscosity end

export ViscosityREFPROP

function scaling_model(param::FrameworkParams{ViscosityREFPROP}, s, x=[1.]) where T
    g = (1.0,1.5,2.0,2.5)
    return generic_powerseries_scaling_model(param, s, x, g)
end

property_CE(prop::ViscosityREFPROP, T, Mw, σ, ε) = viscosity_CE(T, Mw, σ, ε, Ω_22_newfeld(T*kB/ε))

function generic_powerseries_scaling_model(param::FrameworkParams, s, x, g)
    ##α = a1,a2,a3,a4,ξ
    #ln(η + 1) = a1*(s/ξ)^g1 + a2*(s/ξ)^g2 + a3*(s/ξ)^g3 + a4*(s/ξ)^g4
    α = param.α
    res_plus = zero(Base.promote_eltype(α,s,x,g))
    @assert length(g) + 1 <= size(α,1)
    @assert length(x) == size(α,2)
    for i in 1:length(g)
        ni = zero(Base.promote_eltype(α,x))
        for j in 1:length(x)
            ni += x[j]*α[i,j]/(α[end,j]^g[i])
        end
        res_plus += ni*s^g[i]
    end
    return res_plus
end

function scaling(param::FrameworkParams{ViscosityREFPROP}, eos, η⁺p1, T, ϱ, s, z=[1.]; inv=true)
    k = !inv ? 1 : -1
    # ideal viscosity
    ηres = η⁺p1 - 1
    if length(z) == 1
        _1 = one(eltype(z))
        η₀ = _1*property_CE(transport_property(param), T, param.base.Mw[1], param.σ[1], param.ε[1])
    else
        η₀_all = property_CE.(transport_property(param), T, param.base.Mw, param.σ, param.ε)
        η₀ = mix_CE(param.base, η₀_all, z)
    end

    return scaling_property(param, eos, ηres, η₀, T, ϱ, s, z; inv)
end

function scaling_property(param::FrameworkParams{ViscosityREFPROP},eos, ηres, η₀, T, ϱ, s, z; inv = true)
    η₀ + plus_scaling(param.base, ηres, T, ϱ, s, z; inv=inv)
end

calculate_Ymin(::ViscosityREFPROP) = false