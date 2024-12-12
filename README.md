[![Dev][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] [![Build Status][build-img]][build-url]

<p align="center">
  <img width="150px" src="docs/src/assets/logo.svg">
</p>

# EntropyScaling.jl

Transport property modeling based on entropy scaling, equations of state (EOS) and more.

This package provides methods to model 
- the viscosity,
- the thermal conductivity, and
- diffusion coefficients.

For the EOS calculations, additional packages need to be imported. Alternatively, custom EOS
functions can be defined. Implementations of EOS models are *not* included in this package.

The documentation of the package can be found [here][docs-stable-url].

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

## Example

```julia
using EntropyScaling, Clapeyron

# Load sample data
(T_exp,ϱ_exp,η_exp) = EntropyScaling.load_sample_data()
data = ViscosityData(T_exp, [], ϱ_exp, η_exp, :unknown)

# Create EOS model
eos_model = PCSAFT("butane")

# Create entropy scaling model (fit of parameters)
model = FrameworkModel(eos_model, [data])

# Calculation of the viscostiy at state
p = 0.1e6                                       # Pa
T = 300.                                        # K
η = viscosity(model, p, T)
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://se-schmitt.github.io/EntropyScaling.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://se-schmitt.github.io/EntropyScaling.jl/dev

[build-img]: https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml/badge.svg?branch=main
[build-url]: https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml?query=branch%3Amain