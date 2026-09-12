## Viscosity of n-butane: fitting the entropy scaling framework to 14 data points
## and predicting isobars across the gas, liquid, and supercritical regions.

using EntropyScaling, Clapeyron, CairoMakie

(T_exp, ϱ_exp, η_exp) = EntropyScaling.load_sample_data()
data = ViscosityData(T_exp, nothing, ϱ_exp, η_exp, :unknown)

model = ESFramework("n-butane", PCSAFT("n-butane"), [data])

fig = Figure(size=(900, 350))
ax1 = Axis(fig[1, 1])
plot_scaling!(ax1, model, data; cprop=:T, markersize=10)

ax = Axis(fig[1, 3]; xlabel="T / K", ylabel="η / (mPa s)", yscale=log10)
T = 130.0:1.0:500.0
for p in (1e5, 1e6, 1e7)
    lines!(ax, T, [viscosity(model, p, Ti) * 1e3 for Ti in T], label="$(round(Int,p/1e5)) bar")
end
axislegend(ax, "p"; position=:rt, framevisible=false)

for (i, a) in enumerate((ax1, ax))
    text!(a, 0, 1; text="($('a'+i-1))", space=:relative, align=(:left, :top), offset=(6, -4), font=:bold)
end

save(joinpath(@__DIR__, "viscosity.png"), fig; px_per_unit=3)
