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

Internally, the molar volume is always calculated first and then the function `EntropyScaling.VT_transport_property(model, V, T, x)` is called (with `[V] = m³`).

```@docs
EntropyScaling.viscosity
EntropyScaling.thermal_conductivity
EntropyScaling.self_diffusion_coefficient
EntropyScaling.MS_diffusion_coefficient
EntropyScaling.fick_diffusion_coefficient
EntropyScaling.inf_diffusion_coefficient
```