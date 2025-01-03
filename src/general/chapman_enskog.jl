export ChapmanEnskogModel

abstract type AbstractChapmanEnskogModel <: AbstractTransportPropertyModel end

# TODO add LJ parameters from Poling to database
"""
    ChapmanEnskogModel <: AbstractTransportPropertyModel

Chapman-Enskog transport properties for the zero-density limit.

## Fields
- `σ::Vector{T}`: Lennard-Jones size parameter (`[σ] = m`)
- `ε::Vector{T}`: Lennard-Jones energy parameter (`[ε] = J`)
- `Mw::Vector{T}`: molar mass (`[Mw] = kg mol⁻¹`)
- `collision::C`: collision integral method (`KimMonroe()` (default) or `Neufeld()`, see [`Ω`](@ref))

## Constructors

- `ChapmanEnskogModel(components; collision_integral=KimMonroe(), ref="", ref_id="")`: database constructor
- `ChapmanEnskogModel(components, σ, ε, Mw; collision_integral=KimMonroe())`: custom parameters constructor

Input arguments can either be single values (pure) or vectors.
The keywords `ref` (short reference) and `ref_id` (DOI or ISBN) enable the specification of the reference.
Currently, parameters from *Poling et al. (2001)* and *Yang et al. (2022)* are in the database.
Mixture properties are calculated according to the models from *Wilke (1950)* (viscosity), *Mason and Saxen (1958)* (thermal conductivity), and *Miller and Carman (1961)* (self-diffusion).

## Example

```julia
using EntropyScaling 

# Construction with custom parameters
σ, ε, Mw = 3.758e-10, 148.6*EntropyScaling.kB, 16.043e-3            # from Poling et al.
model_methane = ChapmanEnskogModel("methane",σ,ε,Mw)

η_mix = viscosity(model_methane, NaN, 300.)
D_mix = self_diffusion_coefficient(model_methane, NaN, 300.)

# Construction from database
model_mix = ChapmanEnskogModel(["butane","methanol"]; ref="Poling et al. (2001)")

η_mix = viscosity(model_mix, NaN, 300., [.5,.5])
D_mix = self_diffusion_coefficient(model_mix, NaN, 300., [.5,.5])  
```

## References
1.  B. E. Poling, J. M. Prausnitz, and J. P. O’Connell: The Properties of Gases and Liquids, 5th, ed. McGraw-Hill, New York (2001).
2.  X. Yang, X. Xiao, M. Thol, M. Richter, and I. H. Bell: Linking Viscosity to Equations of State Using Residual Entropy Scaling Theory, Int. J. Thermophys. 43 (2022) 183, DOI: https://doi.org/10.1007/s10765-022-03096-9.
3.  C. R. Wilke: A Viscosity Equation for Gas Mixtures, The Journal of Chemical Physics 18 (1950) 517–519, DOI: https://doi.org/10.1063/1.1747673.
4.  E. A. Mason and S. C. Saxena: Approximate Formula for the Thermal Conductivity of Gas Mixtures, The Physics of Fluids 1 (1958) 361–369, DOI: https://doi.org/10.1063/1.1724352.
5.  L. Miller and P. C. Carman: Self-Diffusion in Mixtures. Part 4. -- Comparison of Theory and Experiment for Certain Gas Mixtures, Trans. Faraday Soc. 57 (1961) 2143–2150, DOI: https://doi.org/10.1039/TF9615702143.
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
    Mw_db, ε, σ, refs = out 
    if isempty(Mw)
        Mw = Mw_db
    end
    return ChapmanEnskogModel(comps, σ.*1e-10, ε.*kB, Mw, refs, collision_integral)
end

Base.length(model::AbstractChapmanEnskogModel) = length(model.Mw)

# Viscosity
function viscosity(model::ChapmanEnskogModel, p, T, z=Z1)
    return property_CE(Viscosity(), model, T, z)
end
function viscosity(model::ChapmanEnskogModel, T; i=1)
    return 5/16 * √(model.Mw[i]*kB*T/NA/π) / (model.σ[i]^2*Ω(Viscosity(),model,T;i=i))
end

function viscosity_CE_plus(model::ChapmanEnskogModel, eos, T; i=1)
    if length(eos) == 1
        dBdT = second_virial_coefficient_dT(eos,T)/NA
        B = second_virial_coefficient(eos,T)/NA
    else
        x = zeros(Int64,length(model)); x[i] = 1
        dBdT = second_virial_coefficient_dT(eos,T,x)/NA
        B = second_virial_coefficient(eos,T,x)/NA
    end
    return 5/16/√π / (model.σ[i]^2*Ω(Viscosity(),model,T;i=i)) * (T*dBdT+B)^(2/3)
end

# Thermal conductivity
function thermal_conductivity(model::ChapmanEnskogModel, p, T, z=Z1)
    return property_CE(ThermalConductivity(), model, T, z)
end
function thermal_conductivity(model::ChapmanEnskogModel, T; i=1)
    return 75/64 * kB * √(R*T/model.Mw[i]/π) / (model.σ[i]^2*Ω(ThermalConductivity(),model,T;i=i))
end

function thermal_conductivity_CE_plus(model::ChapmanEnskogModel, eos, T; i=1)
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
function self_diffusion_coefficient(model::ChapmanEnskogModel, p, T, z=Z1)
    return property_CE(SelfDiffusionCoefficient(), model, T, z)
end
function self_diffusion_coefficient(model::ChapmanEnskogModel, T; i=1)
    return 3/8 * √(model.Mw[i]*kB*T/NA/π) / (model.σ[i]^2*Ω(SelfDiffusionCoefficient(),model,T;i=i))
end

function self_diffusion_coefficient_CE_plus(model::ChapmanEnskogModel, eos, T; i=1)
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
function MS_diffusion_coefficient(model::ChapmanEnskogModel, p, T, z)
    return MS_diffusion_coefficient(model::ChapmanEnskogModel, T, z)
end
function MS_diffusion_coefficient(model::ChapmanEnskogModel, T, z=Z1)
    length(model) > 2 && throw(error("Currently only applicable to binary mixtures."))
    return 3/8 * √(model.Mw[1]*kB*T/NA/π) / (model.σ[1]^2*Ω(MaxwellStefanDiffusionCoefficient(),model,T;i=1))
end

function MS_diffusion_coefficient_CE_plus(model::ChapmanEnskogModel, eos, T, z = Z1; i=0)
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
    return MS_diffusion_coefficient_CE(model, T)         # x independent
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

function property_CE_plus(prop::P, model::ChapmanEnskogModel, eos, T, z = Z1; i=0) where
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
    Ω(poperty::AbstractTransportProperty, model::AbstractChapmanEnskogModel, T)

Calculates the collision integral for a given `model` and `property` (`Ω₁₁` for diffusion coefficients and `Ω₂₂` for viscosity/thermal conductivity) at the specified temperature `T`.

Two methods are implemented:
- `KimMonroe()`: *Kim and Monroe (2014)* and
- `Neufeld()`: *Neufeld et al. (1972)*

## References
1.  S. U. Kim and C. W. Monroe: High-Accuracy Calculations of Sixteen Collision Integrals for Lennard-Jones (12-6) Gases and Their Interpolation to Parameterize Neon, Argon, and Krypton, Journal of Computational Physics 273 (2014) 358–373, DOI: https://doi.org/10.1016/j.jcp.2014.05.018.
2.  P. D. Neufeld, A. R. Janzen, and R. A. Aziz: Empirical Equations to Calculate 16 of the Transport Collision Integrals Ω(l,s)* for the Lennard‐Jones (12–6) Potential, The Journal of Chemical Physics 57 (1972) 1100–1102, DOI: https://doi.org/10.1063/1.1678363.
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
    Tc,Pc = crit_pure(eos)
    return correspondence_principle(Tc, Pc)
end

# Viscosity, thermal conductivity: Wilke and Mason and Saxena
struct Wilke <: AbstractTransportPropertyMixing end

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

# Self-diffusion: Miller and Carman
struct MillerCarman <: AbstractTransportPropertyMixing end

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
