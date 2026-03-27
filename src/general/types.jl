export Viscosity, ThermalConductivity, SelfDiffusionCoefficient
export InfDiffusionCoefficient, MaxwellStefanDiffusionCoefficient, FickDiffusionCoefficient

abstract type AbstractTransportPropertyModel end
const ATPM = AbstractTransportPropertyModel
abstract type AbstractEntropyScalingModel <: AbstractTransportPropertyModel end
const AESM = AbstractEntropyScalingModel

abstract type ESFrameworkModel <: AbstractEntropyScalingModel end
abstract type RefpropRESModel  <: AbstractEntropyScalingModel end

abstract type AbstractDiluteGasModel <: AbstractTransportPropertyModel end
abstract type ChapmanEnskogModel <: AbstractDiluteGasModel end

Base.length(model::T) where {T <: AbstractTransportPropertyModel} = length(model.components)
Base.length(model::T) where {T <: AbstractEntropyScalingModel} = length(model.eos)

abstract type AbstractParam end
abstract type AbstractEntropyScalingParams <: AbstractParam end
Base.broadcastable(param::AbstractParam) = Ref(param)

abstract type AbstractTransportProperty end
abstract type AbstractViscosity <: AbstractTransportProperty end
abstract type AbstractThermalConductivity <: AbstractTransportProperty end
Base.broadcastable(prop::AbstractTransportProperty) = Ref(prop)
abstract type DiffusionCoefficient <: AbstractTransportProperty end

struct BaseParam{P} <: AbstractParam
    prop::P
    Mw::CL.SingleParam{Float64}
end

function BaseParam(prop::P, components::Vector{String}, Mw::Vector{Float64};
                   sources=nothing) where {P <: AbstractTransportProperty}
    sp = isnothing(sources) ? CL.SingleParam("Mw", components, Mw) :
                              CL.SingleParam("Mw", components, Mw, fill(false,length(Mw)), sources, nothing)
    return BaseParam(prop, sp)
end

Base.length(base::BaseParam) = length(base.Mw)

abstract type AbstractTransportPropertyData end
const ATPD = AbstractTransportPropertyData

abstract type AbstractTransportPropertyMixing end

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

struct FickDiffusionCoefficient <: DiffusionCoefficient end
name(::FickDiffusionCoefficient) = "Fickian diffusion coefficient"
symbol(::FickDiffusionCoefficient) = :Dᵢⱼ
symbol_name(::FickDiffusionCoefficient) = "D_Fick"

# used for general comparisons
transport_compare_type(P1::AbstractTransportProperty, P2::AbstractTransportProperty) = transport_compare_type(typeof(P1), typeof(P2))
transport_compare_type(P1::Type{T}, P2::Type{T}) where T <: AbstractTransportProperty = true
transport_compare_type(P1::Type{T1}, P2::Type{T2}) where {T1 <: AbstractViscosity, T2 <: AbstractViscosity} = true
transport_compare_type(P1::Type{T1}, P2::Type{T2}) where {T1 <: AbstractThermalConductivity, T2 <: AbstractThermalConductivity} = true
transport_compare_type(P1::Type{T1}, P2::Type{T2}) where {T1 <: AbstractTransportProperty, T2 <: AbstractTransportProperty} = false

struct MSDiffusionMatrix{T} <: AbstractMatrix{T}
    val::Vector{T}
end

Base.size(a::MSDiffusionMatrix) = begin
    N = Int((sqrt(8*length(a.val)+1)-1)/2) + 1
    return (N, N)
end
Base.getindex(a::MSDiffusionMatrix{T}, i, j) where T = begin
    i == j && return T(NaN)
    i < j && return a.val[Int(i/2*(2*size(a,1)-1-i)) + (j-i-1)]
    return a[j,i]
end
Base.setindex!(a::MSDiffusionMatrix{T}, x, i) where T = begin
    setindex!(a.val, x, i)
    return nothing
end
Base.setindex!(a::MSDiffusionMatrix{T}, x, i, j) where T = begin
    if i < j
        _i = Int(i/2*(2*size(a,1)-1-i)) + (j-i-1)
        setindex!(a.val, x, _i)
    elseif i > j
        setindex!(a, x, j, i)
    end
    return nothing
end
Base.zero(::Type{MSDiffusionMatrix}, N) = MSDiffusionMatrix(zeros(sum(1:N-1)))
Base.zero(::Type{MSDiffusionMatrix{T}}, N) where T = MSDiffusionMatrix(zeros(T, sum(1:N-1)))
