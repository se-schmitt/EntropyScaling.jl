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
