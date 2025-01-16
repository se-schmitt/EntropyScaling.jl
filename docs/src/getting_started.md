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
of all required thermodynamic properties, in particular the configurationa entropy $s_{\rm conf}$. 
The EOS calculations are not part of this package.
However, there is an extension to the [`Clapeyron.jl`](https://github.com/ClapeyronThermo/Clapeyron.jl) package,
which provides a large number of different thermodynamic models.
The extension is automatically loaded when loading both packages `EntropyScaling.jl` and `Clapeyron.jl`.
Alternatively, custom thermodynamic models can be used by 'dispatching' the functions defined 
in the [`thermo.jl`](https://github.com/se-schmitt/EntropyScaling.jl/blob/main/src/utils/thermo.jl) file to a custom EOS type.

```julia
using EntropyScaling, Clapeyron

eos_model = PCSAFT("n-butane")
model = FrameworkModel(eos_model,Dict(Viscosity() => [[0.;-14.165;13.97;-2.382;0.501;;]]))
η = viscosity(model, 37.21e6, 323.)
```

Using `EntropyScaling.jl` in combination with [`Clapeyron.jl`](https://github.com/ClapeyronThermo/Clapeyron.jl)
is the recommended way.

## Plots

Plotting functionality is available through [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl) (which must be loaded).
For example, this can be used for checking fitted models against scaled (experimental) reference data (see example below).

```@docs
EntropyScaling.plot
```

#### Example

In the following example, entropy scaling models are fitted to quasi-experimental data (from CoolProp) and compared.

```julia
using EntropyScaling, Clapeyron, CoolProp, Plots

sub = "propane"
Ndat = 200
T, p = rand(200.:500.,Ndat), 10.0.^(rand(Ndat).*4 .+ 4)     # define T-p state points
η = [PropsSI("V","T",T[i],"P",p[i],sub) for i in 1:Ndat]    # calculate reference viscosity data

# Create data and model
ηdat = ViscosityData(T,p,[],η)
model_A = FrameworkModel(PCSAFT(sub), [ηdat])
model_B = FrameworkModel(PCSAFT(sub), [ηdat]; opts=FitOptions(what_fit=Dict(Viscosity() => Bool[0,1,0,1,0])))

# Test plot 
plot(model, ηdat; slims=(0,3), cprop=:T, label="Model A (4 parameters fitted)")
plot!(model_B, ηdat; slims=(0,3), cprop=:T, lc=:blue, label="Model B (2 parameters fitted)")
```
![Plot example](./assets/plot_example.svg)
