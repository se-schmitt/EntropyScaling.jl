export Viscosity, ThermalConductivity, SelfDiffusionCoefficient
export InfDiffusionCoefficient, MaxwellStefanDiffusionCoefficient
export viscosity, thermal_conductivity, self_diffusion_coefficient, MS_diffusion_coefficient

abstract type AbstractEntropyScalingModel end

abstract type AbstractParam end
abstract type AbstractEntropyScalingParams <: AbstractParam end
Base.broadcastable(param::Type{<:AbstractParam}) = Ref(param)
#TODO make parameters broadcastable (see above) and remove unnecessary `Ref(...)`

abstract type AbstractTransportProperty end
abstract type DiffusionCoefficient <: AbstractTransportProperty end

abstract type AbstractTransportPropertyData end

struct Reference
    doi::String
    shortref::String
end
Reference() = Reference("", "NA")

struct EOSInfo <: AbstractParam
    name::String
    ref::Vector{Reference}
end

struct BaseParam{P} <: AbstractParam
    prop::P
    solute_name::Union{Missing,String}
    Mw::Vector{Float64}
    param_ref::Vector{Reference}
    N_data::Int
    T_range::Tuple{Number,Number}
    p_range::Tuple{Number,Number}
end

function BaseParam(prop::P, Mw, dat::D; solute=nothing) where 
                   {P <: AbstractTransportProperty, D <: AbstractTransportPropertyData}
    solute_name = isnothing(solute) ? missing : get_components(solute)[1]
    if prop isa InfDiffusionCoefficient
        Mw = [calc_M_CE([Mw[1],get_Mw(solute)[1]])]
    end
    T_range = (minimum(dat.T), maximum(dat.T))
    p_range = (minimum(dat.p), maximum(dat.p))
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
    viscosity(es_model::EntropyScalingModel, eos, p, T, z=[1.]; phase=:unknown)

Viscosity `η(p,T,x)` (`[η] = Pa s`).
"""
viscosity

"""
    ϱT_viscosity(es_model::EntropyScalingModel, eos, ϱ, T, z=[1.])

Viscosity `η(ϱ,T,x)` (`[η] = Pa s`).
"""
ϱT_viscosity

"""
    thermal_conductivity(es_model::EntropyScalingModel, eos, p, T, z=[1.]; phase=:unknown)

Thermal conductivity `λ(p,T,x)` (`[λ] = W m⁻¹ K⁻¹`).
"""
thermal_conductivity

"""
    ϱT_thermal_conductivity(es_model::EntropyScalingModel, eos, ϱ, T, z=[1.])

Thermal conductivity `λ(ϱ,T,x)` (`[λ] = W m⁻¹ K⁻¹`).
"""
ϱT_thermal_conductivity

"""
    self_diffusion_coefficient(es_model::EntropyScalingModel, eos, p, T, z=[1.]; phase=:unknown)

Self-diffusion coefficient `D(p,T,x)` (`[D] = m² s⁻¹`).
"""
self_diffusion_coefficient

"""
    ϱT_self_diffusion_coefficient(es_model::EntropyScalingModel, eos, ϱ, T, z=[1.])

Self-diffusion coefficient `D(ϱ,T,x)` (`[D] = m² s⁻¹`).
"""
ϱT_self_diffusion_coefficient

"""
    MS_diffusion_coefficient(es_model::EntropyScalingModel, eos, p, T, z; phase=:unknown)

Maxwell-Stefan diffusion coefficient `Ð(p,T,x)` (`[Ð] = m² s⁻¹`).
"""
MS_diffusion_coefficient

"""
    ϱT_MS_diffusion_coefficient(es_model::EntropyScalingModel, eos, ϱ, T, z)
    
Maxwell-Stefan diffusion coefficient `Ð(ϱ,T,x)` (`[Ð] = m² s⁻¹`).
"""
ϱT_MS_diffusion_coefficient