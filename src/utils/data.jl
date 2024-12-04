export TransportPropertyData, ViscosityData, ThermalConductivityData, SelfDiffusionCoefficientData, InfDiffusionCoefficientData
export FitOptions

abstract type AbstractTransportPropertyData end

struct TransportPropertyData{P} <: AbstractTransportPropertyData
    prop::P
    N_dat::Int
    T::Vector{Float64}
    p::Vector{Float64}
    ϱ::Vector{Float64}
    Y::Array{Float64}
    phase::Vector{Symbol}
    ref::Vector{Reference}
end
Base.show(io::IO, data::TransportPropertyData) = print(io,"$(typeof(data))\n    $(data.N_dat) data points.")

function TransportPropertyData(prop, T::Vector, p, ϱ, Y::Vector, phase::Symbol, doi::String="", short::String="")
    N_dat = length(T)
    if isempty(p) && !isempty(ϱ)
        p = ones(N_dat)*NaN
    elseif !isempty(p) && isempty(ϱ)
        ϱ = ones(N_dat)*NaN
    else
        error("Either pressure or density must be provided.")
    end
    [length(k) != N_dat && error("All vectors must have the same length.") for k in [p,ϱ,Y]]

    return TransportPropertyData(prop, length(T), T, p, ϱ, Y, repeat([phase],N_dat), repeat([Reference(doi,short)],N_dat))
end
ViscosityData(args...) = TransportPropertyData(Viscosity(), args...)
ThermalConductivityData(args...) = TransportPropertyData(ThermalConductivity(), args...)
SelfDiffusionCoefficientData(args...) = TransportPropertyData(SelfDiffusionCoefficient(), args...)
InfDiffusionCoefficientData(args...) = TransportPropertyData(InfDiffusionCoefficient(), args...)

function collect_data(datasets::Vector{TPD}, prop::AbstractTransportProperty) where TPD <: TransportPropertyData
    (T, p, ϱ, Y, phase, ref) = (Float64[], Float64[], Float64[], Float64[], Symbol[], Reference[])
    N_dat = 0
    for data in filter_datasets(datasets, prop)
        N_dat += data.N_dat
        push!(T, data.T...)
        push!(p, data.p...)
        push!(ϱ, data.ϱ...)
        push!(Y, data.Y...)
        push!(phase, data.phase...)
        push!(ref, data.ref...)
    end
    return TransportPropertyData(prop, N_dat, T, p, ϱ, Y, phase, ref)
end

filter_datasets(datasets::Vector{TPD}, prop) where TPD <: TransportPropertyData = filter(data -> data.prop == prop, datasets)

struct FitOptions
    what_fit::Dict{T,Vector{Bool}} where T<:AbstractTransportProperty
    FitOptions(;what_fit=Dict{AbstractTransportProperty,Vector{Bool}}()) = new(what_fit)
end