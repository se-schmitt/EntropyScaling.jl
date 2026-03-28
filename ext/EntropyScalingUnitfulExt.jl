module EntropyScalingUnitfulExt

using EntropyScaling
using Unitful
import Unitful: Pa, K, W, m, J, mol, s

const ES = EntropyScaling
Z1 = ES.Z1

const ACEM = ES.ChapmanEnskogModel
const AESM = ES.AbstractEntropyScalingModel
const _tph = Union{Symbol,Vector{Symbol}}

# Density units
Unitful.@derived_dimension __MolarDensity Unitful.𝐍 / Unitful.𝐋^3
Unitful.@derived_dimension __MassDensity Unitful.𝐌 / Unitful.𝐋^3
const __DensityKind = Union{__MolarDensity, __MassDensity}
ustrip_ϱ(ϱ::__MolarDensity, z, Mw) = ustrip(mol/m^3, ϱ)
ustrip_ϱ(ϱ::__MassDensity, z, Mw) = ustrip(mol/m^3, ϱ/ES._dot(Mw, z)u"kg/mol")

# Transport property units
Unitful.@derived_dimension __Viscosity Unitful.𝐌 / Unitful.𝐓 / Unitful.𝐋
Unitful.@derived_dimension __ThermalConductivity Unitful.𝐋 * Unitful.𝐌 / Unitful.𝚯 / Unitful.𝐓^3
Unitful.@derived_dimension __DiffusionCoefficient Unitful.𝐋^2 / Unitful.𝐓

Qeltype(in::Type{Unitful.Quantity{T,D,U}}) where {T,D,U} = T

# Dispatch TransportPropertyData and related
for (fn, dim, unit) in [
        (:ViscosityData, __Viscosity, Pa*s),
        (:ThermalConductivityData, __ThermalConductivity, W/K/m),
        (:SelfDiffusionCoefficientData, __DiffusionCoefficient, m^2/s),
        (:InfDiffusionCoefficientData, __DiffusionCoefficient, m^2/s)
    ]

    @eval begin
        function ES.$fn(T::Vector{<:Unitful.Temperature}, p::Vector{<:Unitful.Pressure}, 
                        Y::Vector{<:$dim}, phase::_tph=:unknown)
            _T, _p, _Y = promote(ustrip.(K,T), ustrip.(Pa,p), ustrip.($unit,Y))
            return ES.$fn(_T, _p, nothing, _Y, phase)
        end
        function ES.$fn(T::Vector{<:Unitful.Temperature}, ϱ::Vector{<:__MolarDensity}, 
                        Y::Vector{<:$dim}, phase::_tph=:unknown)
            _T, _ϱ, _Y = promote(ustrip.(u"K",T), ustrip.(mol/m^3,ϱ), ustrip.($unit,Y))
            return ES.$fn(_T, nothing, _ϱ, _Y, phase)
        end
    end
    if dim != __DiffusionCoefficient
        @eval begin
                function ES.TransportPropertyData(T::Vector{<:Unitful.Temperature}, _pϱ::Vector{PR}, 
                                            Y::Vector{<:$dim}, phase::_tph=:unknown) where 
                                            {PR <: Union{__MolarDensity,Unitful.Pressure}}
                return ES.$fn(T, _pϱ, Y, phase)
            end
        end
    end
end

# Implementation for transport property calculation functions
for (fn,unit) in [
        (:viscosity, Pa*s),
        (:thermal_conductivity, W/K/m),
        (:self_diffusion_coefficient, m^2/s),
        (:MS_diffusion_coefficient, m^2/s)
    ]
    ϱT_fn = Symbol(:ϱT_,fn)
    @eval begin
        # Entropy Scaling models
        function ES.$fn(model::AESM, p::Unitful.Pressure, T::Unitful.Temperature, z=Z1; phase=:unknown, output=$unit)
            _p, _T = ustrip(Pa, p), ustrip(K, T)
            _Y = ES.$fn(model, _p, _T, z; phase)*$unit
            return uconvert(output, _Y)
        end
        function ES.$fn(model::AESM, ϱ::__DensityKind, T::Unitful.Temperature, z=Z1; output=$unit)
            x = z./sum(z)
            _ϱ, _T = ustrip_ϱ(ϱ, x, ES.get_Mw(model.eos)), ustrip(K, T)
            _Y = ES.$ϱT_fn(model, _ϱ, _T, x)*$unit
            return uconvert(output, _Y)
        end

        # Chapman-Enskog models
        function ES.$fn(model::ACEM, p, T::Unitful.Temperature, z=Z1; output=$unit)
            _T = ustrip(K, T)
            _Y = ES.$fn(model, NaN, _T, z)*$unit
            return uconvert(output, _Y)
        end
    end
end

end