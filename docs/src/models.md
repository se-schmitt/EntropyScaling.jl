# Models

```@index
Pages = ["models.md"]
```

## Chapman-Enskog Models

Chapman-Enskog models for transport properties at the zero-density limit based on the kinetic gas theory.

"""@docs
EntropyScaling.ChapmanEnskogModel
"""

## Entropy Scaling Models

Entropy scaling makes use of the fact that transport properies can be scaled such that the
scaled transport property $Y^{\rm s}$ is a univariate function of the configurational (or 
residual) entropy $s_{\rm conf}$, i.e. 

$$Y^{\rm s} = Y^{\rm s}\left(s_{\rm conf}\right).$$

Entropy scaling enables the prediction of transport porperties in all fluid phases.

### Available Models

"""@docs
EntropyScaling.FrameworkModel
EntropyScaling.RefpropRESModel
"""

### Fitting Entropy Scaling Parameters

*[Work in progress]*

Some entropy scaling models allow the fitting of substance-specific parameters to experimental data.
Therefore, a unified interface is provided including the handling of the data and the fit options.

**Fitting Procedure**

1. Loading experimental data
2. Fitting
3. Plotting results and saving the parameters

**Example**