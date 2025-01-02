module UnitfulExt

using EntropyScaling
using Unitful
import Unitful: Pa, K, W, m, J, mol

const ES = EntropyScaling

## Define base Units for dispatch

Unitful.@derived_dimension MolarDensity Unitful.𝐍/Unitful.𝐋^3


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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    ϱ_data::AbstractVector{VR},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}}

    if length(p_data) == 0
        return ES.ViscosityData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            phase,
            doi,
            short
        )
    elseif length(ϱ_data) == 0
        return ES.ViscosityData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            phase,
            doi,
            short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Unitful.Pressure}
    return ES.ViscosityData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        phase,
        doi,
        short
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
        T_data::AbstractVector{T},
        ϱ_data::AbstractVector{VR},
        phase::Symbol,
        doi::String="",
        short::String=""
    ) where {T<:Unitful.Temperature, VR<:MolarDensity}

    return ES.ViscosityData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        phase,
        doi,
        short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    ϱ_data::AbstractVector{VR},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}}

    if length(p_data) == 0
        return ES.ThermalConductivityData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            phase,
            doi,
            short
        )
    elseif length(ϱ_data) == 0
        return ES.ThermalConductivityData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            phase,
            doi,
            short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Unitful.Pressure}

    return ES.ThermalConductivityData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        phase,
        doi,
        short
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
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, VR<:MolarDensity}

    return ES.ThermalConductivityData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        phase,
        doi,
        short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    ϱ_data::AbstractVector{VR},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}}

    if length(p_data) == 0
        return ES.SelfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            phase,
            doi,
            short
        )
    elseif length(ϱ_data) == 0
        return ES.SelfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            phase,
            doi,
            short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Unitful.Pressure}

    return ES.SelfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        phase,
        doi,
        short
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
    T_data::AbstractVector{T},
    ϱ_data::AbstractVector{VR},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, VR<:MolarDensity}

    return ES.SelfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        phase,
        doi,
        short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    ϱ_data::AbstractVector{VR},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Union{Unitful.Pressure,Any},VR<:Union{MolarDensity,Any}}

    if length(p_data) == 0
        return ES.InfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data,
            ϱ_data .|> mol/m^3 .|> ustrip,
            phase,
            doi,
            short
        )
    elseif length(ϱ_data) == 0
        return ES.InfDiffusionCoefficientData(
            T_data .|> K .|> ustrip,
            p_data .|> Pa .|> ustrip,
            ϱ_data,
            phase,
            doi,
            short
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
    T_data::AbstractVector{T},
    p_data::AbstractVector{P},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, P<:Unitful.Pressure}

    return ES.InfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        p_data .|> Pa .|> ustrip,
        [],
        phase,
        doi,
        short
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
    T_data::AbstractVector{T},
    ϱ_data::AbstractVector{VR},
    phase::Symbol,
    doi::String="",
    short::String=""
) where {T<:Unitful.Temperature, VR<:MolarDensity}

    return ES.InfDiffusionCoefficientData(
        T_data .|> K .|> ustrip,
        [],
        ϱ_data .|> mol/m^3 .|> ustrip,
        phase,
        doi,
        short
    )

end



end
