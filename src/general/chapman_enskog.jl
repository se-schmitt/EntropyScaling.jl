export ChapmanEnskogModel, CEPlusModel

abstract type AbstractChapmanEnskogModel <: AbstractTransportPropertyModel end

"""
    ChapmanEnskogModel
"""
struct ChapmanEnskogModel{T,C} <: AbstractChapmanEnskogModel
    σ::Vector{T} 
    ε::Vector{T}
    Mw::Vector{T}       #TODO differentiate between Mw and M_ref
    collision::C 
end
function ChapmanEnskogModel(σ::Vector{T},ε::Vector{T},Mw::Vector{T};collision_integral=KimMonroe()) where {T} 
    return ChapmanEnskogModel(σ,ε,Mw,collision_integral)
end
function ChapmanEnskogModel(σ::Float64,ε::Float64,Mw::Float64;collision_integral=KimMonroe())
    return ChapmanEnskogModel([σ],[ε],[Mw],collision_integral)
end
Base.length(model::AbstractChapmanEnskogModel) = length(Mw)

"""
    CEPlusModel
"""
struct CEPlusModel{T,E,P,C} <: AbstractChapmanEnskogModel
    σ::Vector{T}
    ε::Vector{T}
    Mw::Vector{T}
    eos::E
    collision::C 
end
CEPlusModel(σ,ε,eos,col=KimMonroe()) = CEPlusModel(σ,ε,get_Mw(eos),eos,col)

"""
    viscosity_CE(model::ChapmanEnskogModel, T, z=[1.])

Chapman-Enskog viscosity for the zero-density limit.
"""
function viscosity(model::ChapmanEnskogModel, T; i=1)
    return 5/16 * √(model.Mw[i]*kB*T/NA/π) / (model.σ[i]^2*Ω(model,Viscosity(),T;i=i))
end

"""
    viscosity_plus(model::CEPlusModel, T, z=[1.])

Scaled Chapman-Enskog viscosity for the zero-density limit.
"""
function viscosity_plus(model::CEPlusModel, T; i=1)
    x = zeros(Int64,length(model)); x[i] = 1
    dBdT = second_virial_coefficient_dT(model.eos,T,x)/NA
    B = second_virial_coefficient(model.eos,T,x)/NA
    return 5/16/√π / (model.σ[i]^2*Ω(model,Viscosity(),T;i=i)) * (T*dBdT+B)^(2/3)
end

"""
    thermal_conductivity(model::ChapmanEnskogModel, T, z=[1.])

Chapman-Enskog thermal conductivity for the zero-density limit.
"""
function thermal_conductivity(model::ChapmanEnskogModel, T; i=1)
    return 75/64 * kB * √(R*T/model.Mw[i]/π) / (model.σ[i]^2*Ω(model,ThermalConductivity(),T;i=i))
end

"""
    thermal_conductivity_plus(model:CEPlusModel, T, z=[1.])

Scaled Chapman-Enskog thermal conductivity for the zero-density limit.
"""
function thermal_conductivity_plus(model::CEPlusModel, T; i=1)
    x = zeros(Int64,length(model)); x[i] = 1
    dBdT = second_virial_coefficient_dT(model.eos,T,x)/NA
    B = second_virial_coefficient(model.eos,T,x)/NA
    return 75/64/√π / (model.σ[i]^2*Ω(model,ThermalConductivity(),T;i=i)) * (T*dBdT+B)^(2/3)
end

"""
    self_diffusion_coefficient(model, T, z = [1.0])

Chapman-Enskog diffusion coefficient for the zero-density limit.
"""
function self_diffusion_coefficient(model::ChapmanEnskogModel, T; i=1)
    return 3/8 * √(model.Mw[i]*kB*T/NA/π) / (model.σ[i]^2*Ω(model,DiffusionCoefficient(),T;i=i))
end

"""
    self_diffusion_coefficient_plus(model::CEPlusModel, T, z=[1.])

Scaled Chapman-Enskog self-diffusion coefficient for the zero-density limit.
"""
function self_diffusion_coefficient_plus(model::CEPlusModel, T; i=1)
    x = zeros(Int64,length(model)); x[i] = 1
    dBdT = second_virial_coefficient_dT(model.eos,T,x)/NA
    B = second_virial_coefficient(model.eos,T,x)/NA
    return 3/8/√π / (model.σ[i]^2*Ω(model,DiffusionCoefficient(),T;i=i)) * (T*dBdT+B)^(2/3)
end

"""
    MS_diffusion_coefficient(model::ChapmanEnskogModel, T, z=[1.])

Chapman-Enskog mutual diffusion coefficient for the zero-density limit.
"""
function MS_diffusion_coefficient(model::ChapmanEnskogModel, T, z=[1.])
    length(model) <= 2 && throw(error("Currently only applicable to binary mixtures."))
    return 3/8 * √(model.Mw[1]*kB*T/NA/π) / (model.σ[1]^2*Ω(model,DiffusionCoefficient(),T;i=1))
end

"""
    MS_diffusion_coefficient_plus(model::CEPlusModel, T, z=[1.])

Scaled Chapman-Enskog mutual diffusion coefficient for the zero-density limit.
"""
function MS_diffusion_coefficient_plus(model::CEPlusModel, T, z = Z1)
    length(model) <= 2 && throw(error("Currently only applicable to binary mixtures."))
    dBdT = second_virial_coefficient_dT(model.eos,T,z)/NA
    B = second_virial_coefficient(model.eos,T,z)/NA
    return 3/8/√π / (model.σ[1]^2*Ω(model,DiffusionCoefficient(),T;i=1)) * (T*dBdT+B)^(2/3)
end

# Chapman-Enskog mixture function
function property_CE(model::ChapmanEnskogModel, prop::AbstractTransportProperty, T, z)
    Y₀_all = [property_CE(model, prop, T; i=i) for i in 1:length(model)]
    Y₀ = mix_CE(prop, model, Y₀_all, z)
    return Y₀
end

function property_CE(model::ChapmanEnskogModel, prop::P, T, z=Z1) where 
                    {P <: Union{MaxwellStefanDiffusionCoefficient,InfDiffusionCoefficient}}
    return MS_diffusion_coefficient_CE(model, T)         # x independent
end

# Chapman-Enskog pure functions 
property_CE(model::ChapmanEnskogModel, prop::Viscosity, T; i) = viscosity(model, T; i=i)
property_CE(model::ChapmanEnskogModel, prop::ThermalConductivity, T; i) = thermal_conductivity(model, T; i=i)
property_CE(model::ChapmanEnskogModel, prop::SelfDiffusionCoefficient, T; i) = self_diffusion_coefficient(model, T; i=i)

# Scaled CE mixture function
function property_CE_plus(model::CEPlusModel, prop::AbstractTransportProperty, T, z)
    Y₀⁺_all = [property_CE_plus(model, prop, T; i=i) for i in 1:length(model)]
    Y₀⁺ = mix_CE(prop, model, Y₀⁺_all, z)
    return Y₀⁺
end

function property_CE_plus(model::CEPlusModel, prop::P, T, z = Z1) where
                    {P <: Union{MaxwellStefanDiffusionCoefficient,InfDiffusionCoefficient}}
    return MS_diffusion_coefficient_plus(model, T, z)
end

# Scaled CE pure functions
property_CE_plus(model::CEPlusModel, prop::Viscosity, T; i) = viscosity_plus(model, T; i=i)
property_CE_plus(model::CEPlusModel, prop::ThermalConductivity, T; i) = thermal_conductivity_plus(model, T; i=i)
property_CE_plus(model::CEPlusModel, prop::SelfDiffusionCoefficient, T; i) = self_diffusion_coefficient_plus(model, T; i=i)

# Collision integrals
abstract type AbstractCollisionIntegralMethod end
struct KimMonroe <: AbstractCollisionIntegralMethod end
struct Neufeld <: AbstractCollisionIntegralMethod end

Ω(prop::Union{Viscosity,ThermalConductivity},method::KimMonroe,T_red) = Ω_22(T_red)
Ω(prop::DiffusionCoefficient,method::KimMonroe,T_red) = Ω_11(T_red)
Ω(prop::Union{Viscosity,ThermalConductivity},method::Neufeld,T_red) = Ω_22_neufeld(T_red)
Ω(prop::DiffusionCoefficient,method::Neufeld,T_red) = Ω_11_neufeld(T_red)

Ω(model::AbstractChapmanEnskogModel, prop, T; i) = Ω(prop,model.collision,T*kB/model.ε[i])

"""
    Ω_22(T_red)

Collision integrals from [Kim and Monroe (2014)](https://www.doi.org/10.1016/j.jcp.2014.05.018).
"""
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

"""
    Ω_11(T_red)

Collision integrals from [Kim and Monroe (2014)](https://www.doi.org/10.1016/j.jcp.2014.05.018).
"""
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

"""
    Ω_22_neufeld(T_red)

Collision integrals from [Neufeld (1972)](https://doi.org/10.1063/1.1678363). The non-polynomial terms are neglected (as REFPROP 10.0 does)
"""
function Ω_22_neufeld(T_red)
    return 1.6145*(T_red)^(-0.14874) + 0.52487*exp(-0.77320*T_red) + 2.16178*exp(-2.43787*T_red)
end

"""
    correspondence_principle(Tc, pc)
    correspondence_principle(eos)

Calculate the LJ parameters from the critical temperature and pressure.
"""
function correspondence_principle(Tc, pc)
    Tc_LJ = 1.321
    pc_LJ = 0.129

    ε = kB*Tc/Tc_LJ
    σ = cbrt(pc_LJ/pc*ε)
    return σ, ε
end

function correspondence_principle(eos)
    Tc,Pc = crit_pure(eos)
    return correspondence_principle(Tc, Pc)
end

# Mixing rule from Wilke (1950) [DOI: 10.1063/1.1747673] and Mason and Saxena (1958) [DOI: 10.1063/1.1724352]
struct Wilke <: AbstractTransportPropertyMixing end
"""
    mix_CE(::Wilke,param::BaseParam, Y, x)
    mix_CE(param::BaseParam{Viscosity}, Y, x)
    mix_CE(param::BaseParam{ThermalConductivity}, Y, x)

Mixing rule for Chapman-Enskog transport properties by Wilke (1950) for Viscosity and by 
Mason and Saxena (1958) for ThermalConductivity.

## References
(1) Wilke, C. R. A Viscosity Equation for Gas Mixtures. The Journal of Chemical Physics 
1950, 18 (4), 517–519. https://doi.org/10.1063/1.1747673.
(2) Mason, E. A.; Saxena, S. C. Approximate Formula for the Thermal Conductivity of Gas 
Mixtures. The Physics of Fluids 1958, 1 (5), 361–369. https://doi.org/10.1063/1.1724352.
"""
function mix_CE(::Wilke,model::AbstractChapmanEnskogModel, Y, x)
    Y₀_mix = zero(Base.promote_eltype(Y,x))
    enum_M = enumerate(model.Mw)
    for (i,Mi) in enum_M
        xΦ = zero(Y₀_mix)
        for (j,Mj) in enum_M
            xΦ += x[j] * (1+√(Y[i]/Y[j])*√√(Mj/Mi))^2 / √(8*(1+Mi/Mj))
        end
        Y₀_mix += x[i]*Y[i]/xΦ
    end
    return Y₀_mix
end

struct MillerCarman <: AbstractTransportPropertyMixing end
"""
    mix_CE(::MillerCarman,param::BaseParam, Y, x)
    mix_CE(param::BaseParam{DiffusionCoefficient}, Y, x)

Mixing rule for Chapman-Enskog diffusion coefficient by Miller and Carman (1961).

## References
(1) Miller, L.; Carman, P. C. Self-Diffusion in Mixtures. Part 4. -- Comparison of Theory 
and Experiment for Certain Gas Mixtures. Trans. Faraday Soc. 1961, 57 (0), 2143–2150. 
    https://doi.org/10.1039/TF9615702143.
"""
function mix_CE(::MillerCarman,model::AbstractChapmanEnskogModel, Y, x)
    return 1.0 / sum(x[i] / Y[i] for i in eachindex(Y))
end

function mix_CE(prop::DiffusionCoefficient, model::AbstractChapmanEnskogModel,Y,x)
    return mix_CE(MillerCarman(),model,Y,x)
end

function mix_CE(prop::Union{Viscosity, ThermalConductivity}, model::AbstractChapmanEnskogModel, Y, x)
    return mix_CE(Wilke(),model,Y,x)
end

calc_M_CE(Mw) = 2.0/sum(inv,Mw)
