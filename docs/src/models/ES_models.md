# Entropy Scaling Models

Entropy scaling makes use of the fact that transport properies can be scaled such that the
scaled transport property $Y^{\rm s}$ is a univariate function of the configurational (or 
residual) entropy $s_{\rm conf}$, i.e. 

$$Y^{\rm s} = Y^{\rm s}\left(s_{\rm conf}\right).$$

Entropy scaling enables the prediction of transport porperties in all fluid phases.

The following entropy scaling models are currently implemented:

- Entropy Scaling Framework [schmitt_entropy_2024](@cite) ([`FrameworkModel`](@ref))
- Refprop Residual Entropy Scaling (RES) Model [yang_linking_2022](@cite) ([`RefpropRESModel`](@ref))

All models are similarly structured with the following fields:

- `components::Vector{String}`: names of the chemical components of the system
- `params::Vector{ModelParams}`: vector of model-specific paramater objects
- `eos`: EOS model

The `ModelParams` are model-specific types containing all required parameters of the model.
They always contain the Chapman-Enskog model (`CE_model`) as well as base parameters (`base`)
which itself contains general parameters like the transport  property or the molar mass.

All models share the contructor method `Model(eos, params::Dict{P})`, where params is a dict 
containing the parameters with the respective transport property as key, e.g., 
`Dict(Viscosity() => [a_η, b_η, c_η], ThermalConductivity() => [a_λ, b_λ, c_λ])`.
Here, `a`, `b`, and `c` are the parameters (note that `a_η`, `b_η`, ... are vectors or matrices themselfs).
Lists of the parameters are given below in the repective 'Parameters' sections.
Additional model-specific constructors are also given below.

## Framework Model
```@docs
EntropyScaling.FrameworkModel
```

## Refprop RES Model
```@docs
EntropyScaling.RefpropRESModel
```

## Fitting Utilities

Some entropy scaling models allow the fitting of substance-specific parameters to experimental data.
Therefore, a unified interface is provided including the handling of the data and the fit options.

**Fitting Procedure**

1. Loading experimental data and defining `TransportPropertyData`
2. Fitting (included in the model construction)
3. Plotting results and saving the parameters

```@docs
EntropyScaling.TransportPropertyData
EntropyScaling.FitOptions
```
