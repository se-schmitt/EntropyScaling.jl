# Transport Properties

After constructing a model, transport properties can be calculated (if respective parameters are available).
The call for all this is always of the form:

`transport_property(model, p, T, x=[1.]; phase=:unkwown)`

where
- `model` is a `AbstractTransportPropertyModel`,
- `p` is the pressure (`[p] = Pa`),
- `T` is the temperature (`[T] = K`),
- `x` is the mole fraction (`[x] = mol mol⁻¹`), and
- `phase` is the desired phase (`liquid` or `gas`) used in the volume solver.

Internally, the density is always calculated first and then the function `EntropyScaling.ϱT_transport_property(model, ϱ, T, x)` is called (with `[ϱ] = mol m⁻³`).

```@docs
EntropyScaling.viscosity
EntropyScaling.thermal_conductivity
EntropyScaling.self_diffusion_coefficient
EntropyScaling.MS_diffusion_coefficient
```