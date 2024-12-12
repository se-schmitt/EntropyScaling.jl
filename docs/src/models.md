# Models

```@index
Pages = ["models.md"]
```

## Chapman-Enskog Models

Chapman-Enskog models for transport properties at the zero-density limit based on the kinetic gas theory.

"""@docs
EntropyScaling.viscosity_CE
EntropyScaling.thermal_conductivity_CE
EntropyScaling.diffusion_coefficient_CE
EntropyScaling.viscosity_CE_plus
EntropyScaling.thermal_conductivity_CE_plus
EntropyScaling.diffusion_coefficient_CE_plus
EntropyScaling.Ω_11
EntropyScaling.Ω_22
"""

## Entropy Scaling Models

Entropy scaling makes use of the fact that transport properies can be scaled such that the
scaled transport property $Y^{\rm s}$ is a univariate function of the configurational (or 
residual) entropy $s_{\rm conf}$, i.e. 

$$Y^{\rm s} = Y^{\rm s}\left(s_{\rm conf}\right).$$

Entropy scaling enables the prediction of transport porperties in all fluid phases based on 
few experimental data.

The following entropy scaling methods are implemented:
- [Entropy Scaling Framework](https://doi.org/10.1016/j.molliq.2023.123811)

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


"""@docs
EntropyScaling.FrameworkModel
"""