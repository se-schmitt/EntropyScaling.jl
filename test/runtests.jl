using Test
using TestItemRunner

@testsnippet Plotting begin
    using EntropyScaling, Clapeyron

    # Build a tiny model and small synthetic dataset (viscosity)
    model = ESFramework("n-butane", PCSAFT("n-butane"))

    T = [300.0, 325.0, 350.0]
    p = [1e5, 5e5, 1e6]
    η = [1e-3, 1.5e-3, 2.0e-3]
    ηdat = ViscosityData(T, p, nothing, η)

    # Common keyword sets to exercise
    kw_base = (
        marker=:diamond,
        markersize=7,
        markercolor=:red,
        linecolor=:black,
        linewidth=2,
        linestyle=:dash,
        label="model",
    )
end

@run_package_tests