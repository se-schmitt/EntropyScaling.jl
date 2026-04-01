export ChapmanEnskog

"""
    ChapmanEnskog <: ChapmanEnskogModel
    ChapmanEnskog(components; userlocations=String[], collision_integral=KimMonroe())

Chapman-Enskog transport properties for the zero-density limit.

# Parameters
- `σ::SingleParam{T}`: Lennard-Jones size parameter (`[σ] = Å`)
- `ε::SingleParam{T}`: Lennard-Jones energy parameter (`[ε] = K`)
- `Mw::SingleParam{T}`: molar mass (`[Mw] = g mol⁻¹`)
- `collision::C`: collision integral method (`KimMonroe()` (default) or `Neufeld()`, see [`Ω`](@ref))

Currently, parameters from [poling_properties_2001](@citet) and [yang_linking_2022](@citet) are in the database.
Mixture properties are calculated according to the models from [wilke_viscosity_1950](@citet) (viscosity), [mason_approximate_1958](@citet) (thermal conductivity), and [miller_self-diffusion_1961](@citet) (self-diffusion).

# Example

```julia
using EntropyScaling

# Construction with custom parameters
model_methane = ChapmanEnskog("methane"; userlocations=(;sigma=3.758, epsilon=148.6, Mw=16.043)

η = viscosity(model_methane, NaN, 300.)
D = self_diffusion_coefficient(model_methane, NaN, 300.)

# Construction from database
model_mix = ChapmanEnskog(["butane","methanol"])

η_mix = viscosity(model_mix, NaN, 300., [.5,.5])
D_mix = self_diffusion_coefficient(model_mix, NaN, 300., [.5,.5])
```
"""
struct ChapmanEnskog{T,C} <: ChapmanEnskogModel
    components::Vector{String}
    sigma::CL.SingleParam{T}
    epsilon::CL.SingleParam{T}
    Mw::CL.SingleParam{T}
    collision::C
end

function ChapmanEnskog(components; userlocations=String[], collision_integral=KimMonroe())
    _components = CL.format_components(components)

    params = CL.getparams(_components, [get_db_path(ChapmanEnskogModel, nothing)]; userlocations)

    σ = params["sigma"]
    σ.values .*= 1e-10
    ε = params["epsilon"]
    ε.values .*= kB
    Mw = params["Mw"]
    
    return ChapmanEnskog(_components, σ, ε, Mw, collision_integral)
end

Base.length(model::ChapmanEnskogModel) = length(model.components)

# Viscosity
function viscosity(model::ChapmanEnskogModel, p, T, z=Z1)
    return property_CE(Viscosity(), model, T, z)
end
function viscosity(model::ChapmanEnskogModel, T; i=1)
    σ = model.sigma.values
    Mw = model.Mw.values .* 1e-3
    return 5/16 * √(Mw[i]*kB*T/NA/π) / (σ[i]^2*Ω(Viscosity(),model,T;i=i))
end

function viscosity_CE_plus(model::ChapmanEnskogModel, eos, T; i=1)
    σ = model.sigma.values
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 5/16/√π / (σ[i]^2*Ω(Viscosity(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Thermal conductivity
function thermal_conductivity(model::ChapmanEnskogModel, p, T, z=Z1)
    return property_CE(ThermalConductivity(), model, T, z)
end
function thermal_conductivity(model::ChapmanEnskogModel, T; i=1)
    σ = model.sigma.values
    Mw = model.Mw.values .* 1e-3
    return 75/64 * kB * √(R*T/Mw[i]/π) / (σ[i]^2*Ω(ThermalConductivity(),model,T;i=i))
end

function thermal_conductivity_CE_plus(model::ChapmanEnskogModel, eos, T; i=1)
    σ = model.sigma.values
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 75/64/√π / (σ[i]^2*Ω(ThermalConductivity(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Self-diffusion coefficient
function self_diffusion_coefficient(model::ChapmanEnskogModel, p, T, z=Z1)
    return property_CE(SelfDiffusionCoefficient(), model, T, z)
end
function self_diffusion_coefficient(model::ChapmanEnskogModel, T; i=1)
    σ = model.sigma.values
    Mw = model.Mw.values .* 1e-3
    return 3/8 * √(Mw[i]*kB*T/NA/π) / (σ[i]^2*Ω(SelfDiffusionCoefficient(),model,T;i=i))
end

function self_diffusion_coefficient_CE_plus(model::ChapmanEnskogModel, eos, T; i=1)
    σ = model.sigma.values
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 3/8/√π / (σ[i]^2*Ω(SelfDiffusionCoefficient(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Maxwell-Stefan diffusion coefficient
function MS_diffusion_coefficient(model::ChapmanEnskogModel, p, T, z)
    return MS_diffusion_coefficient(model, T, z)
end
function MS_diffusion_coefficient(model::ChapmanEnskogModel, T, z=Z1)
    σ = model.sigma.values
    Mw = model.Mw.values .* 1e-3
    length(model) > 2 && throw(error("Currently only applicable to binary mixtures."))
    return 3/8 * √(Mw[1]*kB*T/NA/π) / (σ[1]^2*Ω(MaxwellStefanDiffusionCoefficient(),model,T;i=1))
end

function MS_diffusion_coefficient_CE_plus(model::ChapmanEnskogModel, eos, T, z = Z1; i=0)
    length(model) > 2 && throw(error("Currently only applicable to binary mixtures."))
    if i != 0
        z = zeros(Int64,length(eos))
        z[i] = 1
    end

    σ = model.sigma.values
    dBdT = second_virial_coefficient_dT(eos,T,z)/NA
    B = second_virial_coefficient(eos,T,z)/NA
    return 3/8/√π / (σ[1]^2*Ω(MaxwellStefanDiffusionCoefficient(),model,T;i=1)) * (T*dBdT+B)^(2/3)
end

MS_diffusion_coefficient_CE(model, T) = MS_diffusion_coefficient(model, T)

# Chapman-Enskog mixture function
function property_CE(prop::AbstractTransportProperty, model::ChapmanEnskogModel, T, z)
    N_c = length(model)
    if N_c == 1
        return property_CE(prop, model, T; i=1)
    else
        Y₀_all = [property_CE(prop, model, T; i=i) for i in 1:N_c]
        return mix_CE(prop, model, Y₀_all, z)
    end
end

function property_CE(prop::P, model::ChapmanEnskogModel, T, z=Z1) where
                    {P <: Union{MaxwellStefanDiffusionCoefficient,InfDiffusionCoefficient}}
    return MS_diffusion_coefficient_CE(model, T)
end

# Chapman-Enskog pure functions
property_CE(prop::Viscosity, model::ChapmanEnskogModel, T; i) = viscosity(model, T; i=i)
property_CE(prop::ThermalConductivity, model::ChapmanEnskogModel, T; i) = thermal_conductivity(model, T; i=i)
property_CE(prop::SelfDiffusionCoefficient, model::ChapmanEnskogModel, T; i) = self_diffusion_coefficient(model, T; i=i)

# Scaled CE mixture function
function property_CE_plus(prop::AbstractTransportProperty, model::ChapmanEnskogModel, eos, T, z)
    N_c = length(model)
    if N_c == 1
        return property_CE_plus(prop, model, eos, T; i=1)
    else
        Y₀⁺_all = [property_CE_plus(prop, model, eos, T; i=i) for i in 1:N_c]
        return mix_CE(prop, model, Y₀⁺_all, z)
    end
end

function property_CE_plus(prop::P, model::ChapmanEnskogModel, eos, T, z=Z1; i=0) where
                    {P <: Union{MaxwellStefanDiffusionCoefficient,InfDiffusionCoefficient}}
    return MS_diffusion_coefficient_CE_plus(model, eos, T, z; i=i)
end

# Scaled CE pure functions
property_CE_plus(prop::Viscosity, model::ChapmanEnskogModel, eos, T; i) = viscosity_CE_plus(model, eos, T; i=i)
property_CE_plus(prop::ThermalConductivity, model::ChapmanEnskogModel, eos, T; i) = thermal_conductivity_CE_plus(model, eos, T; i=i)
property_CE_plus(prop::SelfDiffusionCoefficient, model::ChapmanEnskogModel, eos, T; i) = self_diffusion_coefficient_CE_plus(model, eos, T; i=i)

# Collision integrals
abstract type AbstractCollisionIntegralMethod end
struct KimMonroe <: AbstractCollisionIntegralMethod end
struct Neufeld <: AbstractCollisionIntegralMethod end

Ω(prop::Union{Viscosity,ThermalConductivity},method::KimMonroe,T_red) = Ω_22(T_red)
Ω(prop::DiffusionCoefficient,method::KimMonroe,T_red) = Ω_11(T_red)
Ω(prop::Union{Viscosity,ThermalConductivity},method::Neufeld,T_red) = Ω_22_neufeld(T_red)
Ω(prop::DiffusionCoefficient,method::Neufeld,T_red) = Ω_11_neufeld(T_red)

"""
    Ω(poperty::AbstractTransportProperty, model::ChapmanEnskogModel, T)

Calculates the collision integral for a given `model` and `property` (`Ω₁₁` for diffusion coefficients and `Ω₂₂` for viscosity/thermal conductivity) at the specified temperature `T`.

Two methods are implemented:
- `KimMonroe()` [kim_high-accuracy_2014](@cite)
- `Neufeld()` [neufeld_empirical_1972](@cite)
"""
Ω(prop::AbstractTransportProperty, model::ChapmanEnskogModel, T; i) = Ω(prop,model.collision,T*kB/model.epsilon[i])

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

function Ω_11_neufeld(T_red)
    return 1.0548*(T_red)^(-0.15504) + 0.54762*exp(-0.77320*T_red) + 2.16178*exp(-2.43787*T_red)
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

# Viscosity, thermal conductivity: Wilke and Mason and Saxena
struct Wilke <: AbstractTransportPropertyMixing end
struct MasonSaxena <: AbstractTransportPropertyMixing end

function mix_CE(::Union{Wilke,MasonSaxena}, model::AbstractDiluteGasModel, Y, x; YΦ=Y)
    Y₀_mix = zero(Base.promote_eltype(Y,x))
    Mw = model.Mw.values .* 1e-3
    enum_M = enumerate(Mw)
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

function mix_CE(::MillerCarman,model::AbstractDiluteGasModel, Y, x)
    return 1.0 / sum(x[i] / Y[i] for i in eachindex(Y))
end

function mix_CE(prop::DiffusionCoefficient, model::AbstractDiluteGasModel,Y,x)
    return mix_CE(MillerCarman(),model,Y,x)
end

function mix_CE(prop::Viscosity, model::AbstractDiluteGasModel, Y, x)
    return mix_CE(Wilke(),model,Y,x)
end

function mix_CE(prop::ThermalConductivity, model::AbstractDiluteGasModel, Y, x)
    return mix_CE(MasonSaxena(),model,Y,x)
end

calc_M_CE(Mw) = 2.0/sum(inv,Mw)
