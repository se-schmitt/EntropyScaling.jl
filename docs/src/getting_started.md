# Getting Started with EntropyScaling.jl

## Installation

The package can be installed by:
```julia
Pkg> add EntropyScaling
```
Package mode can reached by typing `]` in REPL.
Then, the module can be loaded by
```julia
using EntropyScaling
```

## EOS Calculations

The calculation of transport properties through entropy scaling is mostly based on 
fundamental EOS (defined in the Helmholtz energy $a$) as they allow the consistent calculation
of all required thermodynamic properties, in particular the configurational entropy $s_{\rm conf}$. 
The EOS calculations are not part of this package.
[`Clapeyron.jl`](https://github.com/ClapeyronThermo/Clapeyron.jl) is used as thermodynamic backend for the EOS calculations providing a large number of different thermodynamic models.

```julia
using EntropyScaling, Clapeyron

eos_model = PCSAFT("n-butane")
model = ESFramework(eos_model, Dict(Viscosity() => [[0.;-14.165;13.97;-2.382;0.501;;]]))
η = viscosity(model, 37.21e6, 323.)
```

## Units

`EntropyScaling.jl` can be used in combination with [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl).
This enables both to obtain unitful properties directly from the models and to use data with associated units for fitting models.

1. Unitful calculation of [transport properties](@ref "Transport Properties"). This adds the methods:
    - `viscosity(model, p::Pressure, T::Temperature, z=[1.]; phase=:unknown, output=default)`
    - `viscosity(model, ϱ::Density, T::Temperature, z=[1.]; output=default)`
    where `::Property` is meta code for unitful values, e.g. `300u"K"` for `T::Temperature`.
    The density can either be mass or molar density.
    The `output` keyword defines the unit of the calculated transport property (see [here](@ref "Transport Properties") for the default values).
2. Defining data for fitting entropy scaling models using units. 
   This supports constructing `TransportPropertyData` using units, 
   e.g. `TransportPropertyData(T::Vector{Temperature}, p::Vector{Pressure}, η::Vector{Viscosity})`.
   The property-specific constructors (`ViscosityData`, ...) also support units.

In the following, both cases are demonstrated:

```julia
using EntropyScaling, Unitful, Clapeyron, CoolProp

# Calculate unitful transport properties
model_ref = RefpropRES("methane")
thermal_conductivity(model_ref, 1u"atm", 80u"°F")
# 0.019529950158013242 W K^-1 m^-1

# Assign units to data for fitting
(_T_exp,_ϱ_exp,_η_exp) = EntropyScaling.load_sample_data();     # Load sample data
T_exp = _T_exp .* 1u"K"
ϱ_exp = _ϱ_exp .* 1u"mol/m^3"
η_exp = _η_exp .* 1u"Pa*s"

data = TransportPropertyData(T_exp, ϱ_exp, η_exp)

eos_model = PCSAFT("butane")
model = ESFramework(eos_model, [data])                          # Fit model parameters

η = viscosity(model, 1u"bar", 26.85u"°C", phase=:liquid, output = u"cP")
# 0.16058971694885213 cP
```

### Lennard-Jones units

If the equation of state is a Clapeyron potential model (`Thol2016LJ()`, `Johnson1993LJ()`, `KolafaNezbedaLJ()`, `PeTSLJ()`, `Pohl2023Mie6(ν)`, `FeANNMie(ν)`, ...), all calculations are done in Lennard-Jones units (σ = ε = m = k_B = 1).
Inputs are the reduced pressure `p*`, temperature `T*` and density `ϱ*` (or volume `v* = 1/ϱ*`), and the transport properties are returned in reduced form, `η* = η σ²/√(mε)`, `λ* = λ σ² k_B⁻¹ √(m/ε)` and `D* = D σ⁻¹ √(m/ε)`.
The Chapman-Enskog part uses σ* = ε* = 1 by default (custom values can be passed via `σε_ce` in reduced units).

```julia
using EntropyScaling, Clapeyron

eos = Thol2016LJ()
model = ESFramework("LJ", eos; userlocations=Dict(Viscosity() => (; α0=0., α1=1.5, α2=0.3, α3=0., αln=0.)))

viscosity(model, 0.5, 2.0; phase=:liquid)           # η*(p* = 0.5, T* = 2.0)
EntropyScaling.VT_viscosity(model, 1/0.6, 2.0)      # η*(ϱ* = 0.6, T* = 2.0)
viscosity(ChapmanEnskog(eos), NaN, 2.0)             # zero-density limit η₀*(T* = 2.0)
```

## Visualization

Plotting functionality is available through extensions to [`Makie.jl`](https://docs.makie.org/stable/) and [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl).
In particular, "entropy scaling" plots, i.e. the scaled transport property $Y^s$ as function of the scaled entropy $s^s$, can easily be created as described below.

```@docs
EntropyScaling.plot_scaling
```
