# Examples

Scripts reproducing the figures of the JOSS paper.

| Script | Figure | Content |
|---|---|---|
| [`viscosity.jl`](viscosity.jl) | `viscosity.png` | Entropy scaling framework fitted to 14 viscosity data points of *n*-butane; entropy scaling plot and predicted isobars. |
| [`diffusion.jl`](diffusion.jl) | `diffusion.png` | Self-diffusion, Maxwell-Stefan, and Fick diffusion coefficients of *n*-hexane + *n*-dodecane. |

Run from this directory:

```
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
julia --project=. viscosity.jl
julia --project=. diffusion.jl
```
