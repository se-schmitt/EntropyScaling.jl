# Models

```@index
Pages = ["models.md"]
```

## Chapman-Enskog Models

Chapman-Enskog model the for transport properties at the zero-density limit based on the kinetic gas theory.

```math
\begin{aligned}
\eta_{\varrho \rightarrow 0}              &= \frac{5}{16} \sqrt{\frac{M k_{\rm B} T}{\pi N_{\rm A}}} \frac{1}{\sigma^2 \Omega^{(2,2)}} \\
\lambda_{\varrho \rightarrow 0}	          &= \frac{75}{64} k_{\rm B} \sqrt{\frac{R T}{M \pi}} \frac{1}{\sigma^2 \Omega^{(2,2)}}\\
D_{\varrho \rightarrow 0} \varrho^{\rm m} &= \frac{3}{8} \sqrt{\frac{M k_{\rm B} T}{\pi N_{\rm A}}} \frac{1}{\sigma^2 \Omega^{(1,1)}}
\end{aligned}
```

```@docs
EntropyScaling.ChapmanEnskogModel
EntropyScaling.Ω
```

## Entropy Scaling Models

Entropy scaling makes use of the fact that transport properies can be scaled such that the
scaled transport property $Y^{\rm s}$ is a univariate function of the configurational (or 
residual) entropy $s_{\rm conf}$, i.e. 

$$Y^{\rm s} = Y^{\rm s}\left(s_{\rm conf}\right).$$

Entropy scaling enables the prediction of transport porperties in all fluid phases.

### Models

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

```@docs
EntropyScaling.FrameworkModel
EntropyScaling.RefpropRESModel
```

### Fitting Utilities

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