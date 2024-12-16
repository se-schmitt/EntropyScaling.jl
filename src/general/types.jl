export Viscosity, ThermalConductivity, SelfDiffusionCoefficient
export InfDiffusionCoefficient, MaxwellStefanDiffusionCoefficient
export viscosity, thermal_conductivity, self_diffusion_coefficient, MS_diffusion_coefficient

abstract type AbstractTransportPropertyModel end
abstract type AbstractEntropyScalingModel <: AbstractTransportPropertyModel end
Base.length(model::T) where {T<:AbstractEntropyScalingModel} = length(model.eos)

abstract type AbstractParam end
abstract type AbstractEntropyScalingParams <: AbstractParam end
Base.broadcastable(param::AbstractParam) = Ref(param)

abstract type AbstractTransportProperty end
abstract type AbstractViscosity <: AbstractTransportProperty end
abstract type AbstractThermalConductivity <: AbstractTransportProperty end
Base.broadcastable(prop::AbstractTransportProperty) = Ref(prop)
abstract type DiffusionCoefficient <: AbstractTransportProperty end

abstract type AbstractTransportPropertyData end

abstract type AbstractTransportPropertyMixing end

struct Reference
    doi::String
    shortref::String
end
Reference() = Reference("", "NA")

struct BaseParam{P} <: AbstractParam
    prop::P
    solute_name::Union{Missing,String}
    Mw::Vector{Float64}
    param_ref::Vector{Reference}
    N_data::Int
    T_range::Tuple{Float64,Float64}
    p_range::Tuple{Float64,Float64}
end

function BaseParam(prop::P, Mw, dat::D; solute=missing) where
                   {P <: AbstractTransportProperty, D <: AbstractTransportPropertyData}
    solute_name = isnothing(solute) ? missing : get_components(solute)[1]
    if prop isa InfDiffusionCoefficient
        Mw = [calc_M_CE([Mw[1],get_Mw(solute)[1]])]
    end

    T_range = extrema(dat.T)
    p_range = extrema(dat.p)
    return BaseParam(prop, solute_name, Mw, [Reference()], dat.N_dat, T_range, p_range)
end

function BaseParam(prop::P, Mw, ref=[Reference()], N_dat=0, T_range=(NaN,NaN),
                   p_range=(NaN,NaN); solute_name=missing) where
                   {P <: AbstractTransportProperty}


    Mw isa Number && (Mw = [Mw])
    if prop isa InfDiffusionCoefficient
        Mw = [calc_M_CE(Mw) for _ in 1:length(Mw)]
    end
    return BaseParam(prop, solute_name, convert(Vector{Float64},Mw), ref, N_dat, T_range, p_range)
end

Base.length(base::BaseParam) = length(base.Mw)

struct Viscosity <: AbstractViscosity end
name(::AbstractViscosity) = "viscosity"
symbol(::AbstractViscosity) = :η
symbol_name(::AbstractViscosity) = "eta"

struct ThermalConductivity <: AbstractThermalConductivity end
name(::AbstractThermalConductivity) = "thermal conductivity"
symbol(::AbstractThermalConductivity) = :λ
symbol_name(::AbstractThermalConductivity) = "lambda"

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

#used for general comparisons
transport_compare_type(P1::AbstractTransportProperty,P2::AbstractTransportProperty) = transport_compare_type(typeof(P1),typeof(P2))
transport_compare_type(P1::Type{T},P2::Type{T}) where T <: AbstractTransportProperty = true
#get viscosities, thermal conductivities, TODO: do the same structure with diffusion coeffficients?
transport_compare_type(P1::Type{T1},P2::Type{T2}) where {T1 <: AbstractViscosity,T2 <: AbstractViscosity} = true
transport_compare_type(P1::Type{T1},P2::Type{T2}) where {T1 <: AbstractThermalConductivity,T2 <: AbstractThermalConductivity} = true
#fallback
transport_compare_type(P1::Type{T1},P2::Type{T2}) where {T1 <: AbstractTransportProperty,T2 <: AbstractTransportProperty} = false


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

transport_property(x::AbstractTransportProperty) = x
transport_property(x::BaseParam) = x.prop