[![Dev][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] [![Build Status][build-img]][build-url]

<p align="center">
  <img width="150px" src="docs/src/assets/logo.svg">
</p>

# EntropyScaling.jl

Transport property modeling based on entropy scaling, equations of state (EOS) and more.

This package provides models for 
- the viscosity,
- the thermal conductivity, and
- diffusion coefficients.

Currently, the following models are available:
- Chapman-Enskog models for gases
- Entropy Scaling framework

The documentation of the package can be found [here][docs-stable-url].

For the EOS calculations, additional packages need to be imported. Alternatively, custom EOS
functions can be defined. Implementations of EOS models are *not* included in this package.

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

**Chapman-Enskog viscosity of methane**

```julia
using EntropyScaling

# Parameters from Poling et al. (2001)
σ = 3.758e-10                   # size parameter ([σ] = m)
ε = 148.6*EntropyScaling.kB     # energy parameter ([ε] = J)
Mw = 16.0425e-3                 # molar mass ([Mw] = kg/mol)

# Calculate gas viscosity of methane at 300 K
T = 300.
η = viscosity_CE(T, Mw, σ, ε)
```

**Entropy scaling framework in combination with [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)**

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