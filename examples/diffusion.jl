## Diffusion coefficients of the binary mixture n-hexane + n-dodecane at 298.15 K and 1 bar.
## Self-diffusion parameters of the pure components and the binary interaction parameters
## yield the Maxwell-Stefan and (with the thermodynamic factor from the EOS) the Fick
## diffusion coefficient of the mixture.

using EntropyScaling, Clapeyron, CairoMakie

comps = ["hexane", "dodecane"]
model = ESFramework(comps, PCSAFT(comps); userlocations=Dict(
    DiffusionCoefficient() => (; α2=[-2.5414 -4.0610; -2.5463 -4.6610],
                                 α3=[-1.9186 -3.4267; -1.8070 -3.2021]),
))

p, T = 1e5, 298.15
x = 0.0:0.01:1.0
z = [[xi, 1 - xi] for xi in x]

D = [self_diffusion_coefficient(model, p, T, zi) for zi in z]
Ð = [MS_diffusion_coefficient(model, p, T, zi)[1, 2] for zi in z]
D_Fick = [only(fick_diffusion_coefficient(model, p, T, zi)) for zi in z]

fig = Figure(size=(500, 380), fontsize=15)
ax = Axis(fig[1, 1]; xlabel=rich("x", subscript("hexane"), " / (mol mol⁻¹)"),
          ylabel="D / (10⁻⁹ m² s⁻¹)")
lines!(ax, x, first.(D) .* 1e9; linestyle=:dot, label=rich("D", subscript("self,hexane")))
lines!(ax, x, last.(D) .* 1e9; linestyle=:dot, label=rich("D", subscript("self,dodecane")))
lines!(ax, x, Ð .* 1e9; linestyle=:dash, label="Ð (Maxwell-Stefan)")
lines!(ax, x, D_Fick .* 1e9; linewidth=3, label="D (Fick)")
axislegend(ax; position=:lt, framevisible=false)

save(joinpath(@__DIR__, "diffusion.png"), fig; px_per_unit=3)
