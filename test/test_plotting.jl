@testitem "Plotting (Makie.jl)" setup=[Plotting] begin
    using CairoMakie

    # With data, plain markers
    fig = plot_scaling(model, ηdat; kw_base...)
    @test fig isa CairoMakie.Figure

    # With data and color by temperature property and colormap controls
    fig2 = plot_scaling(model, ηdat; kw_base..., cprop=:T, colormap=:viridis, colorrange=(minimum(T), maximum(T)))
    @test fig2 isa CairoMakie.Figure

    # Plot into given axis via plot!
    ax = CairoMakie.Axis(fig2[2,1])
    fig3 = plot_scaling!(ax, model, ηdat; kw_base...)
    @test fig3 isa CairoMakie.Figure

    # No data: provide slims and prop
    fig4 = plot_scaling(model, nothing; slims=(0.0, 3.0), prop=Viscosity(), kw_base...)
    @test fig4 isa CairoMakie.Figure

    # Standalone plot!(model, data) using current axis
    CairoMakie.current_figure() === nothing && (CairoMakie.Figure();)
    fig5 = plot_scaling!(model, ηdat; kw_base...)
    @test fig5 isa CairoMakie.Figure
end

@testitem "Plotting (Plots.jl)" setup=[Plotting] begin
    using Plots

    # With data
    p1 = plot_scaling(model, ηdat; kw_base...)
    @test p1 isa Plots.Plot

    # With data and color by pressure property and colormap controls
    p2 = plot_scaling(model, ηdat; kw_base..., cprop=:p, colormap=:viridis, colorrange=(minimum(p), maximum(p)))
    @test p2 isa Plots.Plot

    # Plot into existing plot via plot!
    plot_scaling!(p2, model, ηdat; kw_base...)
    @test p2 isa Plots.Plot

    # No data: provide slims and prop
    p3 = plot_scaling(model, nothing; slims=(0.0, 3.0), prop=Viscosity(), kw_base...)
    @test p3 isa Plots.Plot

    # Standalone plot!
    p4 = plot_scaling!(model, ηdat; kw_base...)
    @test p4 isa Plots.Plot
end
