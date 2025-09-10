using Test
using EntropyScaling
using Clapeyron

@testset "Plotting" begin
    # Build a tiny model and small synthetic dataset (viscosity)
    eos = PCSAFT("n-butane")
    params = Dict(
        Viscosity() => [[0.0; -14.165; 13.97; -2.382; 0.501;;]],
        ThermalConductivity() => [[3.962; 98.222; -82.974; 20.079; 1.073;;]],
    )
    model = FrameworkModel(eos, params)

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

    @testset "MakieExt" begin
        import CairoMakie

        # With data, plain markers
        fig = CairoMakie.plot(model, ηdat; kw_base...)
        @test fig isa CairoMakie.Figure

        # With data and color by temperature property and colormap controls
        fig2 = CairoMakie.plot(model, ηdat; kw_base..., cprop=:T, colormap=:viridis, colorrange=(minimum(T), maximum(T)))
        @test fig2 isa CairoMakie.Figure

        # Plot into given axis via plot!
        ax = CairoMakie.Axis(fig2[2,1])
        fig3 = CairoMakie.plot!(ax, model, ηdat; kw_base...)
        @test fig3 isa CairoMakie.Figure

        # No data: provide slims and prop
        fig4 = CairoMakie.plot(model, nothing; slims=(0.0, 3.0), prop=Viscosity(), kw_base...)
        @test fig4 isa CairoMakie.Figure

        # Standalone plot!(model, data) using current axis
        CairoMakie.current_figure() === nothing && (CairoMakie.Figure();)
        fig5 = CairoMakie.plot!(model, ηdat; kw_base...)
        @test fig5 isa CairoMakie.Figure
    end

    @testset "PlotsExt" begin
        import Plots

        # With data
        p1 = Plots.plot(model, ηdat; kw_base...)
        @test p1 isa Plots.Plot

        # With data and color by pressure property and colormap controls
        p2 = Plots.plot(model, ηdat; kw_base..., cprop=:p, colormap=:viridis, colorrange=(minimum(p), maximum(p)))
        @test p2 isa Plots.Plot

        # Plot into existing plot via plot!
        Plots.plot!(p2, model, ηdat; kw_base...)
        @test p2 isa Plots.Plot

        # No data: provide slims and prop
        p3 = Plots.plot(model, nothing; slims=(0.0, 3.0), prop=Viscosity(), kw_base...)
        @test p3 isa Plots.Plot

        # Standalone plot!
        p4 = Plots.plot!(model, ηdat; kw_base...)
        @test p4 isa Plots.Plot
    end
end
