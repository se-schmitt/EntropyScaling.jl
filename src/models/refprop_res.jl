export RefpropRESModel, RefpropRESParams

struct RefpropRESParams{P,T} <: AbstractEntropyScalingParams
    n::Matrix{T}
    ξ::Vector{T}
    g::Vector{T}
    crit::Dict{Symbol,Vector{T}}
    CE_model::ChapmanEnskogModel
    base::BaseParam{P}
end

function RefpropRESParams(prop::AbstractTransportProperty, eos, n::Matrix{T}, ξ::Vector{T},
    g::Vector{T}, σ::Vector{T}, ε::Vector{T}, 
    crit::Dict{Symbol,Vector{T}}=Dict{Symbol,Vector{T}}()) where T
    
    Mw = convert(typeof(ξ),get_Mw(eos))
    CE_model = ChapmanEnskogModel(repeat([""],length(Mw)),σ,ε,Mw,collision_integral=KimMonroe())
    base = BaseParam(prop, Mw)
    return RefpropRESParams(n,ξ,g,crit,CE_model,base)
end

"""
    RefpropRESModel{T} <: AbstractEntropyScalingModel

Entropy scaling model based on Refprop EOS [yang_linking_2022,yang_entropy_2021](@cite).

A database provides ready-to-use models for the viscosity of several fluids.
The model can be used in combination with [`Clapeyron.jl`](https://github.com/ClapeyronThermo/Clapeyron.jl) and [`Coolprop.jl`](https://github.com/CoolProp/CoolProp.jl) (see examples).

# Parameters 

- `n::Matrix{T}`: component-specific *or* global (group) parameters
- `ξ::Vector{T}`: component-specific scaling parameter in case global parameters are used (`ξ = 1` for individual fits)
- `σ::Vector{T}`: LJ size parameter for the Chapman-Enskog model 
- `ε::Vector{T}`: LJ energy parameter for the Chapman-Enskog model
- `crit::Dict{Symbol,Vector}`: parameters for critical contribution of thermal conductivity (keys: `:φ0`, `:Γ`, `:qD`, and `:Tref`)

# Constructors

- `RefpropRESModel(eos, params::Dict{P}; kw...)`: Default constructor (see above).
- `RefpropRESModel(eos, components; kw...)`: Creates a ES model using the parameters provided in the database (recommended). 
    `RefpropRESModel(components; kw...)` creates the EOS model on-the-fly (only works if `Clapeyron.jl` and `Coolprop.jl` are loaded).

!!! info
    The default CoolProp EOS is used here which does not necessarily match the choice of the original papers. 
    This might lead to slight deviations to the values in the original papers (especially for the thermal conductivity).

**Keywords**

- `ηref = "Yang et al. (2022)"`: viscosity model (`"Yang et al. (2022)"` [yang_linking_2022](@cite) or `"Martinek et al. (2025)"` [martinek_entropy_2025](@cite)).

# Example

```julia 
using EntropyScaling, Clapeyron, CoolProp

model_pure = RefpropRESModel("R134a")
η_pure = viscosity(model_pure, 1e5, 300.; phase=:liquid)

model_mix = RefpropRESModel(["decane","butane"])
η_mix = viscosity(model_mix, 1e5, 300., [.5,.5])
```
"""
struct RefpropRESModel{E,P} <: AbstractEntropyScalingModel
    components::Vector{<:AbstractString}
    params::P
    eos::E
end

@modelmethods RefpropRESModel RefpropRESParams

function RefpropRESModel(eos, components::Vector{<:AbstractString}; ηref=nothing)
    _ηref = isnothing(ηref) ? "Yang et al. (2022)" : ηref

    params = RefpropRESParams[]
    for prop in [Viscosity(),ThermalConductivity()]
        propref = prop == Viscosity() ? _ηref : ""
        g = (prop == Viscosity() && ηref == "Martinek et al. (2025)") ? [1.8, 2.4, 2.8] : 
            [1.0,1.5,2.0,2.5]
        out = load_params(RefpropRESModel, prop, components; ref=propref)
        if !ismissing(out)
            if prop == ThermalConductivity()
                ξ, n1, n2, n3, n4, φ0, Γ, qD, Tref, refs = out
                crit = Dict(:φ0 => φ0, :Γ => Γ, :qD => qD, :Tref => Tref)
            else
                ξ, n1, n2, n3, n4, refs = out
                crit = Dict{Symbol,Vector{eltype(ξ)}}()
            end
            CE_model = ChapmanEnskogModel(components; Mw=get_Mw(eos), ref_id="10.1007/s10765-022-03096-9")
            base = BaseParam(prop, get_Mw(eos), refs)
            push!(params, RefpropRESParams(permutedims(hcat(n1,n2,n3,n4)),ξ,g,crit,CE_model,base))
        end
    end
    isempty(params) ? throw(MissingException("No parameters found for system [$(join(components,", "))]")) : nothing
    return RefpropRESModel(components, Tuple(params), eos)
end

function scaling_model(param::RefpropRESParams, s, x=z1)
    return exp(powerseries_scaling_model(param, s, x, param.g)) - 1.
end

function scaling_model(param::RefpropRESParams{P}, s, x=z1) where P <: ThermalConductivity
    return powerseries_scaling_model(param, s, x, param.g)
end

function powerseries_scaling_model(param, s, x, g)
    n, ξ = param.n, param.ξ
    Y⁺ = zero(Base.promote_eltype(n,s,x,g))
    @assert length(g) <= size(n,1)
    @assert length(x) == size(n,2)
    for i in eachindex(g)
        ni = zero(Base.promote_eltype(n,x))
        for j in eachindex(x)
            ni += x[j]*n[i,j]/(ξ[j]^g[i])
        end
        Y⁺ += ni*s^g[i]
    end
    return Y⁺ # Y⁺
end

function scaling(param::RefpropRESParams, eos, Y, T, ϱ, s, z=Z1; inv=true, η=nothing)
    k = !inv ? 1 : -1
    tp = transport_property(param)

    if tp == ThermalConductivity()
        each_ind = eachindex(z)
        λ₀ = [thermal_conductivity(param.CE_model, T; i) for i in each_ind]
        η₀ = [viscosity(param.CE_model, T; i) for i in each_ind]
        pure_eos = split_model(eos)
        cₚ₀ = isobaric_heat_capacity.(pure_eos, 1e-10, T)
        λ_int = thermal_conductivity_internal.(η₀, cₚ₀, get_Mw(eos))
        Δλ_c = thermal_conductivity_critical(param.crit, eos, ϱ, T, η, z)
        Y₀ = mix_CE(MasonSaxena(), param.CE_model, λ₀ .+ λ_int, z; YΦ=λ₀) + Δλ_c
    else
        Y₀ = property_CE(tp, param.CE_model, T, z)
    end

    if inv
        return plus_scaling(param.base, Y, T, ϱ, s, z; inv=true) + Y₀
    else
        return plus_scaling(param.base, Y - Y₀, ϱ, s, z; inv=false)
    end
end

function ϱT_thermal_conductivity(model::RefpropRESModel, ϱ, T, z::AbstractVector=Z1)
    param = model[ThermalConductivity()]
    s = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)
    λˢ = scaling_model(param, sˢ, z)
    η = ϱT_viscosity(model, ϱ, T, z)
    return scaling(param, model.eos, λˢ, T, ϱ, s, z; inv=true, η)
end

function thermal_conductivity_internal(η₀, cₚ, M)
    f_int = 1.32
    return f_int*η₀/M * (cₚ - 5/2*R)
end

function thermal_conductivity_critical(crit_param, eos, ϱ, T, η, z)
    # Constants
    RD = 1.02
    ν_γ = 0.63 / 1.239

    # Parameters 
    φ0, Γ, qD, Tref = [_dot(crit_param[k],z) for k in (:φ0, :Γ, :qD, :Tref)]
    
    # EOS
    ∂ϱ∂p = ForwardDiff.derivative(xp -> molar_density(eos, xp, T, z; ϱ0=ϱ), pressure(eos, ϱ, T, z))
    ∂ϱ∂p_Tref = ForwardDiff.derivative(xp -> molar_density(eos, xp, Tref, z; ϱ0=ϱ), pressure(eos, ϱ, Tref, z))
    Δ∂ϱ∂p = ∂ϱ∂p - Tref/T*∂ϱ∂p_Tref
    if Δ∂ϱ∂p < 0
        return 0.0
    else
        cₚ = isobaric_heat_capacity(eos, ϱ, T, z)
        cᵥ = isochoric_heat_capacity(eos, ϱ, T, z)
        _, pc, ϱc = crit_mix(eos, z)
    
        cᵥ_cₚ = cᵥ/cₚ
        ϱr = ϱ/ϱc
        φ = φ0 * (pc/ϱc*ϱr/Γ * Δ∂ϱ∂p)^(ν_γ)
        qDφ = qD*φ
        ΔΩ = 2/π * ((1-cᵥ_cₚ)*atan(qDφ) + cᵥ_cₚ*qDφ - 1 + exp(-qDφ/(1+qDφ/3*(qDφ/ϱr)^2)))
        
        return ϱ*cₚ*RD*kB*T/(6π*η*φ) * ΔΩ
    end
end