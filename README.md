[![Dev][docs-stable-img]][docs-stable-url] [![Build Status][build-img]][build-url]

# EntropyScaling.jl

Transport property modeling based on entropy scaling and equations of state (EOS).

This package provides methods to model 
- the viscosity,
- the thermal conductivity, and
- diffusion coefficients

in a physically sound way. For the EOS calculations, additional packages need to be imported.
Alternatively, custom EOS functions can be defined. Implementations of EOS models are *not*
included in this package.

Entropy scaling makes use of the fact that transport properies can be scaled such that the
scaled transport property $Y^{\rm s}$ is a univariate function of the configurational (or 
residual) entropy $s_{\rm conf}$, i.e. 
$$Y^{\rm s} = Y^{\rm s}\left(s_{\rm conf}\right).$$

Entropy scaling enables the prediction of transport porperties in all fluid phases based on 
few experimental data.

The following entropy scaling methods are implemented:
- [Entropy Scaling Framework](https://doi.org/10.1016/j.molliq.2023.123811)
- more to come...

All methods are based on empirical parameters fitted to experimental data of the respective
transport property. If no parameters are given for a specific substance, the workflow for 
most methods is the following:
1. **Creating the entropy scaling model**: Fitting of empirical parameters to experimental 
   data. 
2. **Calculating transport properties**: Evaluating the entropy scaling model at any fluid
   state point.
In general, parameters are not transferable between different EOS models, i.e. for they 
should only be used (in step 2) in combination with the EOS which was also used for fitting 
the parameters (in step 1).

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

[build-img]: https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml?query=branch%3Amain
[build-url]: https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml/badge.svg?branch=main