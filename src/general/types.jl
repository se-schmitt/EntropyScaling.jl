export Viscosity, ThermalConductivity, SelfDiffusionCoefficient
export InfDiffusionCoefficient, MaxwellStefanDiffusionCoefficient
export viscosity, thermal_conductivity, self_diffusion_coefficient, MS_diffusion_coefficient

abstract type AbstractEntropyScalingModel end
Base.length(model::T) where {T<:AbstractEntropyScalingModel} = length(model.eos)

abstract type AbstractParam end
abstract type AbstractEntropyScalingParams <: AbstractParam end
Base.broadcastable(param::T) where {T<:AbstractParam} = Ref(param)

abstract type AbstractTransportProperty end
Base.broadcastable(prop::T) where {T<:AbstractTransportProperty} = Ref(prop)
abstract type DiffusionCoefficient <: AbstractTransportProperty end

abstract type AbstractTransportPropertyData end

struct Reference
    doi::String
    shortref::String
end
Reference() = Reference("", "NA")

struct BaseParam{P,T} <: AbstractParam
    prop::P
    solute_name::Union{Missing,String}
    Mw::Vector{Float64}
    param_ref::Vector{Reference}
    N_data::Int
    T_range::Tuple{T,T}
    p_range::Tuple{T,T}
end

function BaseParam(prop::P, Mw, dat::D; solute=nothing) where 
                   {P <: AbstractTransportProperty, D <: AbstractTransportPropertyData}
    solute_name = isnothing(solute) ? missing : get_components(solute)[1]
    if prop isa InfDiffusionCoefficient
        Mw = [calc_M_CE([Mw[1],get_Mw(solute)[1]])]
    end
    TT = Base.promote_eltype(dat.T,dat.p)

    T_range = TT.(extrema(dat.T))
    p_range = TT.(extrema(dat.p))
    return BaseParam(prop, solute_name, Mw, [Reference()], dat.N_dat, T_range, p_range)
end

function BaseParam(prop::P, Mw, ref=[Reference()], N_dat=0, T_range=(NaN,NaN), 
                   p_range=(NaN,NaN); solute_name=missing) where 
                   P <: AbstractTransportProperty
    if prop isa InfDiffusionCoefficient
        Mw = [calc_M_CE(Mw) for _ in 1:length(Mw)]
    end
    return BaseParam(prop, solute_name, Mw, ref, N_dat, T_range, p_range)
end

function BaseParam{P,TT}(prop::P, Mw, ref=[Reference()], N_dat=0, T_range=(NaN,NaN), 
    p_range=(NaN,NaN); solute_name=missing) where
    {TT,P <: AbstractTransportProperty}
    if prop isa InfDiffusionCoefficient
        Mw = [calc_M_CE(Mw) for _ in 1:length(Mw)]
    end
    return BaseParam(prop, solute_name, Mw, ref, N_dat, TT.(T_range), TT.(p_range))
end

function Base.convert(::Type{BaseParam{P,T1}},param::BaseParam{P,T2}) where {P,T1,T2}
    prop = param.prop
    solute_name = param.solute_name
    Mw = param.Mw
    param_ref = param.param_ref
    N_data = param.N_data
    T_range = T1.(param.T_range)
    p_range = T1.(param.p_range)
    return BaseParam(prop,solute_name,Mw,param_ref,N_data,T_range,p_range)
end

# function BaseParam(prop::P, )

#     return BaseParam(prop, Mw, [Reference("fit")], dat.N_dat, T_range, p_range)
# end

struct Viscosity <: AbstractTransportProperty end
name(::Viscosity) = "viscosity"
symbol(::Viscosity) = :η
symbol_name(::Viscosity) = "eta"

struct ThermalConductivity <: AbstractTransportProperty end
name(::ThermalConductivity) = "thermal conductivity"
symbol(::ThermalConductivity) = :λ
symbol_name(::ThermalConductivity) = "lambda"

struct SelfDiffusionCoefficient <: DiffusionCoefficient end
name(::SelfDiffusionCoefficient) = "self-diffusion coefficient"
symbol(::SelfDiffusionCoefficient) = :D
symbol_name(::SelfDiffusionCoefficient) = "D"

struct InfDiffusionCoefficient <: DiffusionCoefficient end
name(::InfDiffusionCoefficient) = "diffusion coefficient at infinite dilution"
symbol(::InfDiffusionCoefficient) = :D∞
symbol_name(::InfDiffusionCoefficient) = "D_inf"

struct MaxwellStefanDiffusionCoefficient <: DiffusionCoefficient end
name(::MaxwellStefanDiffusionCoefficient) = "Maxwell-Stefan diffusion coefficient"
symbol(::MaxwellStefanDiffusionCoefficient) = :Ð
symbol_name(::MaxwellStefanDiffusionCoefficient) = "D_MS"

"""
    viscosity(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Viscosity `η(p,T,x)` (`[η] = Pa s`).
"""
viscosity

"""
    ϱT_viscosity(model::EntropyScalingModel, ϱ, T, z=[1.])

Viscosity `η(ϱ,T,x)` (`[η] = Pa s`).
"""
ϱT_viscosity

"""
    thermal_conductivity(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Thermal conductivity `λ(p,T,x)` (`[λ] = W m⁻¹ K⁻¹`).
"""
thermal_conductivity

"""
    ϱT_thermal_conductivity(model::EntropyScalingModel, ϱ, T, z=[1.])

Thermal conductivity `λ(ϱ,T,x)` (`[λ] = W m⁻¹ K⁻¹`).
"""
ϱT_thermal_conductivity

"""
    self_diffusion_coefficient(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Self-diffusion coefficient `D(p,T,x)` (`[D] = m² s⁻¹`).
"""
self_diffusion_coefficient

"""
    ϱT_self_diffusion_coefficient(model::EntropyScalingModel, ϱ, T, z=[1.])

Self-diffusion coefficient `D(ϱ,T,x)` (`[D] = m² s⁻¹`).
"""
ϱT_self_diffusion_coefficient

"""
    MS_diffusion_coefficient(model::EntropyScalingModel, p, T, z; phase=:unknown)

Maxwell-Stefan diffusion coefficient `Ð(p,T,x)` (`[Ð] = m² s⁻¹`).
"""
MS_diffusion_coefficient

"""
    ϱT_MS_diffusion_coefficient(model::EntropyScalingModel, ϱ, T, z)
    
Maxwell-Stefan diffusion coefficient `Ð(ϱ,T,x)` (`[Ð] = m² s⁻¹`).
"""
ϱT_MS_diffusion_coefficient