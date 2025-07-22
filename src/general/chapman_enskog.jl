export ChapmanEnskogModel

abstract type AbstractChapmanEnskogModel <: AbstractTransportPropertyModel end

"""
    ChapmanEnskogModel <: AbstractTransportPropertyModel

Chapman-Enskog transport properties for the zero-density limit.

# Fields
- `σ::Vector{T}`: Lennard-Jones size parameter (`[σ] = m`)
- `ε::Vector{T}`: Lennard-Jones energy parameter (`[ε] = J`)
- `Mw::Vector{T}`: molar mass (`[Mw] = kg mol⁻¹`)
- `collision::C`: collision integral method (`KimMonroe()` (default) or `Neufeld()`, see [`Ω`](@ref))

# Constructors

- `ChapmanEnskogModel(components; collision_integral=KimMonroe(), ref="", ref_id="")`: database constructor
- `ChapmanEnskogModel(components, σ, ε, Mw; collision_integral=KimMonroe())`: custom parameters constructor

Input arguments can either be single values (pure) or vectors.
The keywords `ref` (short reference) and `ref_id` (DOI or ISBN) enable the specification of the reference.
Currently, parameters from [poling_properties_2001](@citet) and [yang_linking_2022](@citet) are in the database.
Mixture properties are calculated according to the models from [wilke_viscosity_1950](@citet) (dynamic_viscosity), [mason_approximate_1958](@citet) (thermal conductivity), and [miller_self-diffusion_1961](@citet) (self-diffusion).

# Example

```julia
using EntropyScaling 

# Construction with custom parameters
σ, ε, Mw = 3.758e-10, 148.6*EntropyScaling.kB, 16.043e-3            # from Poling et al.
model_methane = ChapmanEnskogModel("methane",σ,ε,Mw)

η_mix = dynamic_viscosity(model_methane, NaN, 300.)
D_mix = self_diffusion_coefficient(model_methane, NaN, 300.)

# Construction from database
model_mix = ChapmanEnskogModel(["butane","methanol"]; ref="Poling et al. (2001)")

η_mix = dynamic_viscosity(model_mix, NaN, 300., [.5,.5])
D_mix = self_diffusion_coefficient(model_mix, NaN, 300., [.5,.5])  
```
"""
struct ChapmanEnskogModel{T,C} <: AbstractChapmanEnskogModel
    components::Vector{String}
    σ::Vector{T} 
    ε::Vector{T}
    Mw::Vector{T}
    ref::Vector{Reference}
    collision::C 
end

function ChapmanEnskogModel(comps::Vector{String}, σ::Vector{T}, ε::Vector{T}, Mw::Vector{T}; collision_integral=KimMonroe()) where {T} 
    return ChapmanEnskogModel(comps, σ, ε, Mw, Reference[], collision_integral)
end

function ChapmanEnskogModel(comps::String, σ::Float64, ε::Float64, Mw::Float64; collision_integral=KimMonroe())
    return ChapmanEnskogModel([comps], [σ], [ε], [Mw], Reference[], collision_integral)
end

ChapmanEnskogModel(comps::String; kwargs...) = ChapmanEnskogModel([comps]; kwargs...)
function ChapmanEnskogModel(comps::Vector{String}; Mw=[], ref="", ref_id="", collision_integral=KimMonroe())
    out = load_params(ChapmanEnskogModel, "", comps; ref, ref_id)
    ismissing(out) ? throw(MissingException("No CE parameters found for system [$(join(comps,", "))]")) : nothing
    Mw_db, ε, σ, refs = out 
    if isempty(Mw)
        Mw = Mw_db
    end
    return ChapmanEnskogModel(comps, σ.*1e-10, ε.*kB, Mw, refs, collision_integral)
end

Base.length(model::AbstractChapmanEnskogModel) = length(model.Mw)

# DynamicViscosity
function dynamic_viscosity(model::AbstractChapmanEnskogModel, p, T, z=Z1)
    return property_CE(DynamicViscosity(), model, T, z)
end
function dynamic_viscosity(model::AbstractChapmanEnskogModel, T; i=1)
    return 5/16 * √(model.Mw[i]*kB*T/NA/π) / (model.σ[i]^2*Ω(DynamicViscosity(),model,T;i=i))
end

function dynamic_viscosity_CE_plus(model::AbstractChapmanEnskogModel, eos, T; i=1)
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 5/16/√π / (model.σ[i]^2*Ω(DynamicViscosity(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Thermal conductivity
function thermal_conductivity(model::AbstractChapmanEnskogModel, p, T, z=Z1)
    return property_CE(ThermalConductivity(), model, T, z)
end
function thermal_conductivity(model::AbstractChapmanEnskogModel, T; i=1)
    return 75/64 * kB * √(R*T/model.Mw[i]/π) / (model.σ[i]^2*Ω(ThermalConductivity(),model,T;i=i))
end

function thermal_conductivity_CE_plus(model::AbstractChapmanEnskogModel, eos, T; i=1)
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 75/64/√π / (model.σ[i]^2*Ω(ThermalConductivity(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Self-diffusion coefficient
function self_diffusion_coefficient(model::AbstractChapmanEnskogModel, p, T, z=Z1)
    return property_CE(SelfDiffusionCoefficient(), model, T, z)
end
function self_diffusion_coefficient(model::AbstractChapmanEnskogModel, T; i=1)
    return 3/8 * √(model.Mw[i]*kB*T/NA/π) / (model.σ[i]^2*Ω(SelfDiffusionCoefficient(),model,T;i=i))
end

function self_diffusion_coefficient_CE_plus(model::AbstractChapmanEnskogModel, eos, T; i=1)
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 3/8/√π / (model.σ[i]^2*Ω(SelfDiffusionCoefficient(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Maxwell-Stefan diffusion coefficient
function MS_diffusion_coefficient(model::AbstractChapmanEnskogModel, p, T, z)
    return MS_diffusion_coefficient(model::AbstractChapmanEnskogModel, T, z)
end
function MS_diffusion_coefficient(model::AbstractChapmanEnskogModel, T, z=Z1)
    length(model) > 2 && throw(error("Currently only applicable to binary mixtures."))
    return 3/8 * √(model.Mw[1]*kB*T/NA/π) / (model.σ[1]^2*Ω(MaxwellStefanDiffusionCoefficient(),model,T;i=1))
end

function MS_diffusion_coefficient_CE_plus(model::AbstractChapmanEnskogModel, eos, T, z = Z1; i=0)
    length(model) > 2 && throw(error("Currently only applicable to binary mixtures."))
    if i != 0
        z = zeros(Int64,length(eos))
        z[i] = 1
    end

    dBdT = second_virial_coefficient_dT(eos,T,z)/NA
    B = second_virial_coefficient(eos,T,z)/NA
    return 3/8/√π / (model.σ[1]^2*Ω(MaxwellStefanDiffusionCoefficient(),model,T;i=1)) * (T*dBdT+B)^(2/3)
end

# Chapman-Enskog mixture function
function property_CE(prop::AbstractTransportProperty, model::AbstractChapmanEnskogModel, T, z)
    N_c = length(model)
    if N_c == 1
        return property_CE(prop, model, T; i=1) 
    else
        Y₀_all = [property_CE(prop, model, T; i=i) for i in 1:N_c]
        return mix_CE(prop, model, Y₀_all, z)
    end
end

function property_CE(prop::P, model::AbstractChapmanEnskogModel, T, z=Z1) where 
                    {P <: Union{MaxwellStefanDiffusionCoefficient,InfDiffusionCoefficient}}
    return MS_diffusion_coefficient_CE(model, T)         # x independent
end

# Chapman-Enskog pure functions 
property_CE(prop::DynamicViscosity, model::AbstractChapmanEnskogModel, T; i) = dynamic_viscosity(model, T; i=i)
property_CE(prop::ThermalConductivity, model::AbstractChapmanEnskogModel, T; i) = thermal_conductivity(model, T; i=i)
property_CE(prop::SelfDiffusionCoefficient, model::AbstractChapmanEnskogModel, T; i) = self_diffusion_coefficient(model, T; i=i)

# Scaled CE mixture function
function property_CE_plus(prop::AbstractTransportProperty, model::AbstractChapmanEnskogModel, eos, T, z)
    N_c = length(model)
    if N_c == 1
        return property_CE_plus(prop, model, eos, T; i=1)
    else
        Y₀⁺_all = [property_CE_plus(prop, model, eos, T; i=i) for i in 1:N_c]
        return mix_CE(prop, model, Y₀⁺_all, z)
    end
end

function property_CE_plus(prop::P, model::AbstractChapmanEnskogModel, eos, T, z = Z1; i=0) where
                    {P <: Union{MaxwellStefanDiffusionCoefficient,InfDiffusionCoefficient}}
    return MS_diffusion_coefficient_CE_plus(model, eos, T, z; i=i)
end

# Scaled CE pure functions
property_CE_plus(prop::DynamicViscosity, model::AbstractChapmanEnskogModel, eos, T; i) = dynamic_viscosity_CE_plus(model, eos, T; i=i)
property_CE_plus(prop::ThermalConductivity, model::AbstractChapmanEnskogModel, eos, T; i) = thermal_conductivity_CE_plus(model, eos, T; i=i)
property_CE_plus(prop::SelfDiffusionCoefficient, model::AbstractChapmanEnskogModel, eos, T; i) = self_diffusion_coefficient_CE_plus(model, eos, T; i=i)

# Collision integrals
abstract type AbstractCollisionIntegralMethod end
struct KimMonroe <: AbstractCollisionIntegralMethod end
struct Neufeld <: AbstractCollisionIntegralMethod end

Ω(prop::Union{DynamicViscosity,ThermalConductivity},method::KimMonroe,T_red) = Ω_22(T_red)
Ω(prop::DiffusionCoefficient,method::KimMonroe,T_red) = Ω_11(T_red)
Ω(prop::Union{DynamicViscosity,ThermalConductivity},method::Neufeld,T_red) = Ω_22_neufeld(T_red)
Ω(prop::DiffusionCoefficient,method::Neufeld,T_red) = Ω_11_neufeld(T_red)

"""
    Ω(property::AbstractTransportProperty, model::AbstractChapmanEnskogModel, T)

Calculates the collision integral for a given `model` and `property` (`Ω₁₁` for diffusion coefficients and `Ω₂₂` for dynamic viscosity/thermal conductivity) at the specified temperature `T`.

Two methods are implemented:
- `KimMonroe()` [kim_high-accuracy_2014](@cite)
- `Neufeld()` [neufeld_empirical_1972](@cite)
"""
Ω(prop::AbstractTransportProperty, model::AbstractChapmanEnskogModel, T; i) = Ω(prop,model.collision,T*kB/model.ε[i])

# Ω_22
function Ω_22(T_red)
    A = -0.92032979
    BCi = ( (   2.3508044,      1.6330213      ),
            (   0.50110649,     -6.9795156e-1  ),
            (   -4.7193769e-1,  1.6096572e-1,  ),
            (   1.5806367e-1,   -2.2109440e-2  ),
            (   -2.6367184e-2,  1.7031434e-3   ),
            (   1.8120118e-3,   -0.56699986e-4 ))
    Ω22 = zero(1.0*T_red) + A
    T_red_i = T_red
    lnT_red = log(T_red)
    lnT_red_i = lnT_red
    for k in 1:6
        BCi1,BCi2 = BCi[k]
        Ω22 += BCi1/T_red_i + BCi2*lnT_red_i
        T_red_i *= T_red
        lnT_red_i *= lnT_red
    end
    return Ω22
end

# Ω_11
function Ω_11(T_red)
    A = -1.1036729
    BCi = ( (   2.6431984,      1.6690746      ),
            (   0.0060432255,   -6.9145890e-1  ),
            (   -1.5158773e-1,  1.5502132e-1,  ),
            (   0.54237938e-1,  -2.0642189e-2  ),
            (   -0.90468682e-2, 1.5402077e-3   ),
            (   0.61742007e-3,  -0.49729535e-4 ))
    Ω11 = zero(1.0*T_red) + A
    T_red_i = T_red
    lnT_red = log(T_red)
    lnT_red_i = lnT_red
    for k in 1:6
        BCi1,BCi2 = BCi[k]
        Ω11 += BCi1/T_red_i + BCi2*lnT_red_i
        T_red_i *= T_red
        lnT_red_i *= lnT_red
    end
    return Ω11
end

# Ω_22_neufeld
function Ω_22_neufeld(T_red)
    return 1.6145*(T_red)^(-0.14874) + 0.52487*exp(-0.77320*T_red) + 2.16178*exp(-2.43787*T_red)
end

# Correspondence principle
function correspondence_principle(Tc, pc)
    Tc_LJ = 1.321
    pc_LJ = 0.129

    ε = kB*Tc/Tc_LJ
    σ = cbrt(pc_LJ/pc*ε)
    return σ, ε
end

function correspondence_principle(eos)
    Tc,pc,_ = crit_pure(eos)
    return correspondence_principle(Tc, pc)
end

# Dynamic viscosity, thermal conductivity: Wilke and Mason and Saxena
struct Wilke <: AbstractTransportPropertyMixing end
struct MasonSaxena <: AbstractTransportPropertyMixing end

function mix_CE(::Union{Wilke,MasonSaxena}, model::AbstractChapmanEnskogModel, Y, x; YΦ=Y)
    Y₀_mix = zero(Base.promote_eltype(Y,x))
    enum_M = enumerate(model.Mw)
    for (i,Mi) in enum_M
        xΦ = zero(Y₀_mix)
        for (j,Mj) in enum_M
            xΦ += x[j] * (1+√(YΦ[i]/YΦ[j])*√√(Mj/Mi))^2 / √(8*(1+Mi/Mj))
        end
        Y₀_mix += x[i]*Y[i]/xΦ
    end
    return Y₀_mix
end

# Self-diffusion: Miller and Carman
struct MillerCarman <: AbstractTransportPropertyMixing end

function mix_CE(::MillerCarman,model::AbstractChapmanEnskogModel, Y, x)
    return 1.0 / sum(x[i] / Y[i] for i in eachindex(Y))
end

function mix_CE(prop::DiffusionCoefficient, model::AbstractChapmanEnskogModel,Y,x)
    return mix_CE(MillerCarman(),model,Y,x)
end

function mix_CE(prop::DynamicViscosity, model::AbstractChapmanEnskogModel, Y, x)
    return mix_CE(Wilke(),model,Y,x)
end

function mix_CE(prop::ThermalConductivity, model::AbstractChapmanEnskogModel, Y, x)
    return mix_CE(MasonSaxena(),model,Y,x)
end

calc_M_CE(Mw) = 2.0/sum(inv,Mw)
