export TransportPropertyData, ViscosityData, ThermalConductivityData
export SelfDiffusionCoefficientData, InfDiffusionCoefficientData
export FitOptions

struct TransportPropertyData{P} <: AbstractTransportPropertyData
    prop::P
    N_dat::Int
    T::Vector{Float64}
    p::Vector{Float64}
    ϱ::Vector{Float64}
    Y::Array{Float64}
    phase::Vector{Symbol}
end
function Base.show(io::IO, data::TransportPropertyData)
    print(io,"$(typeof(data))\n $(data.N_dat) data points.")
    return nothing
end

"""
    TransportPropertyData(prop, T, p, ϱ, η, phase=:unknown)
    ViscosityData(T, p, ϱ, η, phase=:unknown)
    ThermalConductivityData(T, p, ϱ, λ, phase=:unknown)
    SelfDiffusionCoefficientData(T, p, ϱ, D, phase=:unknown)
    InfDiffusionCoefficientData(T, p, ϱ, D, phase=:unknown)

Constructor for `TransportPropertyData`.
Either pressure `p` or density `ϱ` have to be specified.
`phase` can also be a `Vector{Symbol}`.

## Units

- `[T] = K`
- `[p] = Pa`
- `[ϱ] = mol m⁻³`
- `[η] = Pa s`
- `[λ] = W (m K)⁻¹`
- `[D] = m² s⁻¹`
"""
function TransportPropertyData(prop, T::Vector, p, ϱ, Y::Vector,
                               phase::Union{Symbol,Vector{Symbol}}=:unknown)
    
    N_dat = length(T)
    if isempty(p) && !isempty(ϱ)
        p = ones(N_dat)*NaN
    elseif !isempty(p) && isempty(ϱ)
        ϱ = ones(N_dat)*NaN
    else
        error("Either pressure or density must be provided.")
    end

    if phase isa Symbol
        phases = repeat([phase],N_dat)
    else
        phases = phase 
    end

    for k in [p,ϱ,Y,phases]
        length(k) != N_dat && 
            throw(DimensionMismatch("All vectors must have the same length."))
    end

    return TransportPropertyData(prop, length(T), T, p, ϱ, Y, phases)
end
ViscosityData(args...) = TransportPropertyData(Viscosity(), args...)
ThermalConductivityData(args...) = TransportPropertyData(ThermalConductivity(), args...)
SelfDiffusionCoefficientData(args...) = TransportPropertyData(SelfDiffusionCoefficient(), args...)
InfDiffusionCoefficientData(args...) = TransportPropertyData(InfDiffusionCoefficient(), args...)

function collect_data(  datasets::Vector{TPD}, prop::AbstractTransportProperty) where
                        TPD <: TransportPropertyData
    (T, p, ϱ, Y, phase) = (Float64[], Float64[], Float64[], Float64[], Symbol[])
    N_dat = 0
    for data in filter_datasets(datasets, prop)
        N_dat += data.N_dat
        push!(T, data.T...)
        push!(p, data.p...)
        push!(ϱ, data.ϱ...)
        push!(Y, data.Y...)
        push!(phase, data.phase...)
    end
    return TransportPropertyData(prop, N_dat, T, p, ϱ, Y, phase)
end

function filter_datasets(datasets::Vector{TPD}, prop) where TPD <: TransportPropertyData
    return filter(data -> data.prop == prop, datasets)
end

"""
    FitOptions

Struct to control fitting.

## Fields 

- `what_fit::Dict{AbstractTransportProperty,Vector{Bool}}`: specify which parameters to fit

## Example
```
FitOptions(;
    what_fit=Dict(
        ThermalConductivity()=>ones(Bool,5), 
        SelfDiffusionCoefficient()=>Bool[0,0,0,1,1])
)
```
"""
struct FitOptions
    what_fit::Dict{T,Vector{Bool}} where T<:AbstractTransportProperty
    FitOptions(;what_fit=Dict{AbstractTransportProperty,Vector{Bool}}()) = new(what_fit)
end