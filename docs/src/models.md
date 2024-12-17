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

Entropy scaling enables the prediction of transport porperties in all fluid phases based on 
few experimental data.

The following entropy scaling methods are implemented:
- [Entropy Scaling Framework](https://doi.org/10.1016/j.molliq.2023.123811)
- [Refprop RES Model](https://doi.org/10.1007/s10765-022-03096-9)

"""@docs
EntropyScaling.FrameworkModel
EntropyScaling.RefpropRESModel
"""