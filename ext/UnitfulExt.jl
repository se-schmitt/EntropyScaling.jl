module UnitfulExt

using EntropyScaling
using Unitful
import Unitful: Pa, K, W, m, J, mol, s

const ES = EntropyScaling

## Define base Units for dispatch

Unitful.@derived_dimension MolarDensity Unitful.𝐍/Unitful.𝐋^3
Unitful.@derived_dimension ThermalCond Unitful.𝐋*Unitful.𝐌/Unitful.𝚯/Unitful.𝐓^3
Unitful.@derived_dimension Visc Unitful.𝐌/Unitful.𝐓/Unitful.𝐋
Unitful.@derived_dimension DiffusionCoefficient Unitful.𝐋^2/Unitful.𝐓




## Dispatch on the Unit types


"""
    ES.ViscosityData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a ViscosityData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures possible to be empty array.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities possible to be empty array.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.ViscosityData(
    T_data::Vector{T},
    p_data::Vector{P},
    ϱ_data::Vector{VR},
    η_data::Vector{Eta},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},Eta<:Visc,VR<:Union{MolarDensity,Any}}

    if length(p_data) == 0
        return ES.ViscosityData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            η_data .|> Pa*s .|> ustrip,
            phase,
        )
    elseif length(ϱ_data) == 0
        return ES.ViscosityData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            η_data .|> Pa*s .|> ustrip,
            phase,
        )
    else
        error("Either pressure or density must be provided.")
    end
end


"""
    ES.ViscosityData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, phase::Symbol, doi::String="", short::String="")

    Create a ViscosityData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.ViscosityData(
    T_data::Vector{T},
    p_data::Vector{P},
    η_data::Vector{Eta},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Unitful.Pressure, Eta<:Visc}
    return ES.ViscosityData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        η_data .|> Pa*s .|> ustrip,
        phase,
    )
end

"""
    ES.ViscosityData(T_data::AbstractVector{T}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a ViscosityData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.
"""
function ES.ViscosityData(
        T_data::Vector{T},
        ϱ_data::Vector{VR},
        η_data::Vector{Eta},
        phase::Union{Symbol,Vector{Symbol}}=:unknown,
    ) where {T<:Unitful.Temperature, VR<:MolarDensity, Eta<:Visc}

    return ES.ViscosityData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        η_data .|> Pa*s .|> ustrip,
        phase,
    )
end


"""
    ES.ThermalConductivityData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a ThermalConductivityData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures possible to be empty array.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities possible to be empty array.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.ThermalConductivityData(
    T_data::Vector{T},
    p_data::Vector{P},
    ϱ_data::Vector{VR},
    λ_data::Vector{TC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}, TC<:ThermalCond}

    if length(p_data) == 0
        return ES.ThermalConductivityData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            λ_data .|> W/(m*K) .|> ustrip,
            phase,
        )
    elseif length(ϱ_data) == 0
        return ES.ThermalConductivityData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            λ_data .|> W/(m*K) .|> ustrip,
            phase,
        )
    else
        error("Either pressure or density must be provided.")
    end

end

"""
    ES.ThermalConductivityData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, phase::Symbol, doi::String="", short::String="")

    Create a ThermalConductivityData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.ThermalConductivityData(
    T_data::Vector{T},
    p_data::Vector{P},
    λ_data::Vector{TC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Unitful.Pressure, TC<:ThermalCond}

    return ES.ThermalConductivityData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        λ_data .|> W/(m*K) .|> ustrip,
        phase,
    )

end

"""
    ES.ThermalConductivityData(T_data::AbstractVector{T}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a ThermalConductivityData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.ThermalConductivityData(
    T_data::AbstractVector{T},
    ϱ_data::AbstractVector{VR},
    λ_data::AbstractVector{TC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, VR<:MolarDensity, TC<:ThermalCond}

    return ES.ThermalConductivityData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        λ_data .|> W/(m*K) .|> ustrip,
        phase,
    )

end

"""
    ES.SelfDiffusionCoefficientData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a SelfDiffusionCoefficientData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures possible to be empty array.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities possible to be empty array.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.SelfDiffusionCoefficientData(
    T_data::Vector{T},
    p_data::Vector{P},
    ϱ_data::Vector{VR},
    D_data::Vector{DC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}, DC<:DiffusionCoefficient}

    if length(p_data) == 0
        return ES.SelfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            D_data .|> m^2/s .|> ustrip,
            phase,
        )
    elseif length(ϱ_data) == 0
        return ES.SelfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            D_data .|> m^2/s .|> ustrip,
            phase,
        )
    else
        error("Either pressure or density must be provided.")
    end

end

"""
    ES.SelfDiffusionCoefficientData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, phase::Symbol, doi::String="", short::String="")

    Create a SelfDiffusionCoefficientData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.SelfDiffusionCoefficientData(
    T_data::Vector{T},
    p_data::Vector{P},
    D_data::Vector{DC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Unitful.Pressure, DC<:DiffusionCoefficient}

    return ES.SelfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        D_data .|> m^2/s .|> ustrip,
        phase,
    )

end

"""
    ES.SelfDiffusionCoefficientData(T_data::AbstractVector{T}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a SelfDiffusionCoefficientData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.SelfDiffusionCoefficientData(
    T_data::Vector{T},
    ϱ_data::Vector{VR},
    D_data::Vector{DC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, VR<:MolarDensity, DC<:DiffusionCoefficient}

    return ES.SelfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        D_data .|> m^2/s .|> ustrip,
        phase,
    )

end

"""
    ES.InfDiffusionCoefficientData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a InfDiffusionCoefficientData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures possible to be empty array.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities possible to be empty array.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.InfDiffusionCoefficientData(
    T_data::Vector{T},
    p_data::Vector{P},
    ϱ_data::Vector{VR},
    D_data::Vector{DC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}, DC<:DiffusionCoefficient}

    if length(p_data) == 0
        return ES.InfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            D_data .|> m^2/s .|> ustrip,
            phase,
        )
    elseif length(ϱ_data) == 0
        return ES.InfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            D_data .|> m^2/s .|> ustrip,
            phase,
        )
    else
        error("Either pressure or density must be provided.")
    end

end

"""
    ES.InfDiffusionCoefficientData(T_data::AbstractVector{T}, p_data::AbstractVector{P}, phase::Symbol, doi::String="", short::String="")

    Create a InfDiffusionCoefficientData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `p_data::AbstractVector{P}`: Vector of pressures.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.InfDiffusionCoefficientData(
    T_data::Vector{T},
    p_data::Vector{P},
    D_data::Vector{DC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, P<:Unitful.Pressure, DC<:DiffusionCoefficient}

    return ES.InfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        D_data .|> m^2/s .|> ustrip,
        phase,
    )

end

"""
    ES.InfDiffusionCoefficientData(T_data::AbstractVector{T}, ϱ_data::AbstractVector{VR}, phase::Symbol, doi::String="", short::String="")

    Create a InfDiffusionCoefficientData object from the given data.

    # Arguments
    - `T_data::AbstractVector{T}`: Vector of temperatures.
    - `ϱ_data::AbstractVector{VR}`: Vector of molar densities.
    - `phase::Symbol`: Phase of the data.

    # Returns
    - `TransportPropertyData`: TransportPropertyData struct.

"""
function ES.InfDiffusionCoefficientData(
    T_data::Vector{T},
    ϱ_data::Vector{VR},
    D_data::Vector{DC},
    phase::Union{Symbol,Vector{Symbol}}=:unknown,
) where {T<:Unitful.Temperature, VR<:MolarDensity, DC<:DiffusionCoefficient}

    return ES.InfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        D_data .|> m^2/s .|> ustrip,
        phase,
    )

end



end
