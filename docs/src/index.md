# EntropyScaling.jl

Transport property modeling based on entropy scaling and equations of state (EoS).

This package provides methods to model 
- the dynamic viscosity,
- the thermal conductivity, and
- diffusion coefficients

in a physically sound way. For the EoS calculations, additional packages need to be imported.
Alternatively, custom EoS functions can be defined. Implementations of EoS models are *not*
included in this package.

Entropy scaling makes use of the fact that transport properties can be scaled such that the
scaled transport property $Y^{\rm s}$ is a univariate function of the configurational (or 
residual) entropy $s_{\rm conf}$, i.e. 
$$Y^{\rm s} = Y^{\rm s}\left(s_{\rm conf}\right).$$

Entropy scaling enables the prediction of transport properties in all fluid phases based on 
few experimental data.
