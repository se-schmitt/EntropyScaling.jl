@testitem "internals" begin
    using Clapeyron

    @testset "split_model" begin
        @testset "ESFramework" begin
            # Match the working mixture pattern from test_framework.jl: PCSAFT has no
            # benzene/hexane data, so thermal_conductivity is supplied via userlocations.
            # That keeps the (separately-buggy) multi-component diffusion branch dormant.
            comps = ["benzene", "hexane"]
            model = ESFramework(comps, PCSAFT(comps);
                userlocations=Dict(
                    ThermalConductivity() => (;
                        α0=[6.492054, 9.756965],
                        α2=[1.985553, 1.657211],
                        α3=[3.126438, 5.058428],
                    ),
                ),
            )

            pures = Clapeyron.split_model(model)
            @test pures isa AbstractVector
            @test length(pures) == 2
            @test all(m -> m isa typeof(model).name.wrapper, pures)
            @test all(m -> length(m) == 1, pures)
            for (i, m) in pairs(pures)
                @test m.components == [comps[i]]
                @test m.sources == model.sources
                @test length(m.eos) == 1
                @test m.eos.components == [comps[i]]
            end

            # Pure submodel evaluates to the same thermal conductivity as the mixture
            # at pure composition.
            p, T = 0.1e6, 294.7
            for i in 1:2
                z = zeros(2); z[i] = 1.0
                λ_split = thermal_conductivity(pures[i], p, T)
                λ_mix = thermal_conductivity(model, p, T, z)
                @test isfinite(λ_split)
                @test λ_split ≈ λ_mix rtol=1e-10
            end

            # Explicit splitter: keep both components as a single 2-component submodel
            sub = Clapeyron.split_model(model, [[1, 2]])
            @test length(sub) == 1
            @test length(sub[1]) == 2
            @test sub[1].components == comps
        end

        @testset "RefpropRES" begin
            # methane + octane: both have RefpropRES viscosity params AND are in
            # Clapeyron's bundled Empiric database (no CoolProp dep needed).
            comps = ["methane", "octane"]
            model = RefpropRES(comps)
            pures = Clapeyron.split_model(model)
            @test length(pures) == 2
            @test pures[1].components == ["methane"]
            @test pures[2].components == ["octane"]
            @test pures[1].sources == model.sources

            p, T = 1e5, 300.0
            η_split = viscosity(pures[1], p, T; phase=:liquid)
            η_ref = viscosity(RefpropRES("methane"), p, T; phase=:liquid)
            @test η_split ≈ η_ref rtol=1e-10
        end

        @testset "GCES" begin
            model = GCES(["hexane", "benzene"])
            pures = Clapeyron.split_model(model)
            @test length(pures) == 2
            @test pures[1].components == ["hexane"]
            @test pures[2].components == ["benzene"]
            @test length(pures[1].eos) == 1
            @test length(pures[2].eos) == 1

            p, T = 1e6, 300.0
            @test isfinite(viscosity(pures[1], p, T))
            @test isfinite(viscosity(pures[2], p, T))
        end
    end
end
