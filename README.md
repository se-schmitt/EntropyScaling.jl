[![Dev][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] [![Build Status][build-img]][build-url] [![DOI][zenodo-img]][zenodo-url]

<p align="center">
  <img width="150px" src="docs/src/assets/logo.svg">
</p>

# EntropyScaling.jl

Transport property modeling based on entropy scaling and equations of state (EOS).

This package provides models for 
- the viscosity,
- the thermal conductivity, and
- diffusion coefficients.

The documentation of the package can be found [here][docs-dev-url].

For the EOS calculations, additional packages need to be imported. Alternatively, custom EOS
functions might be defined. Implementations of EOS models are *not* included in this package.

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

## Examples

**Chapman-Enskog viscosity (zero-density limit) of methane**

```julia
julia> using EntropyScaling

julia> model = ChapmanEnskog("methane")
ChapmanEnskogModel{methane}
 σ: [3.758] Å
 ε: [148.6] K
 M: [0.01604] kg/m³
 Collision integral: KimMonroe

julia> η = viscosity(model, NaN, 300.)      # gas viscosity of methane at 300 K in Pa s
1.1189976373570321e-5
```

**Fitting a new entropy scaling model**

```julia
julia> using EntropyScaling, Clapeyron

julia> (T_exp,ϱ_exp,η_exp) = EntropyScaling.load_sample_data();    # Load sample data
┌ Info: Experimental data for the viscosity of n-butane.
└       Units: [T] = K, [ϱ] = mol/m³, [η] = Pa·s

julia> data = ViscosityData(T_exp, [], ϱ_exp, η_exp, :unknown)
TransportPropertyData{Viscosity}
    15 data points.

julia> eos_model = PCSAFT("butane")                     # Clapeyron.jl EOS model
PCSAFT{BasicIdeal, Float64} with 1 component:
 "butane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model = FrameworkModel(eos_model, [data])        # Fit model parameters
FrameworkModel with 1 component:
 "butane"
 Available properties: viscosity
 Equation of state: Clapeyron.EoSVectorParam{PCSAFT{BasicIdeal, Float64}}("butane")

julia> η = viscosity(model, 1e5, 300.; phase=:liquid)   # viscostiy at p=1 bar and T=300 K
0.0001605897169488518
```

**Integration with Unitful.jl**

It is possible to work woth units using [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl).
The following example uses the model from above, but calculates the viscosity with units.

```julia
julia> using Unitful

julia> η = viscosity(model, 1u"bar", 26.85u"°C", phase=:liquid, output = u"cP")
0.16058971694885213 cP
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://se-schmitt.github.io/EntropyScaling.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://se-schmitt.github.io/EntropyScaling.jl/dev

[build-img]: https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml/badge.svg?branch=main
[build-url]: https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml?query=branch%3Amain

[zenodo-img]: https://zenodo.org/badge/723050048.svg
[zenodo-url]: https://zenodo.org/badge/latestdoi/723050048
