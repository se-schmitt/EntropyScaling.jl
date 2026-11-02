@testitem "Potential models (LJ units)" begin
    using Clapeyron
    ES = EntropyScaling
    eos = Thol2016LJ()
    coeffs = Dict(
        Viscosity() => (; α0=0., α1=1.5, α2=0.3, α3=0., αln=0.),
        ThermalConductivity() => (; α0=1., α1=1.2, α2=0.2, α3=0., αln=0.),
        SelfDiffusionCoefficient() => (; α0=zeros(1,1), α1=fill(-1.,1,1), α2=fill(-.2,1,1), α3=zeros(1,1), αln=zeros(1,1)))

    @testset "Chapman-Enskog" begin
        ce = ChapmanEnskog(eos)
        for T in (0.8, 2.0, 5.0)
            η = viscosity(ce, NaN, T)
            @test η * ES.Ω(Viscosity(), ce, T; i=1) ≈ 5/16 * sqrt(T/π)
            @test thermal_conductivity(ce, NaN, T) / η ≈ 15/4
        end
    end

    @testset "ESFramework equals SI LJRef in reduced units" begin
        T, ϱ = 2.0, 0.6
        es = ESFramework("LJ", eos; userlocations=coeffs)
        Y = (ES.VT_viscosity(es, 1/ϱ, T), ES.VT_thermal_conductivity(es, 1/ϱ, T),
            ES.VT_self_diffusion_coefficient(es, 1/ϱ, T)[1])
        @test viscosity(es, pressure(eos, 1/ϱ, T), T; phase=:liquid) ≈ Y[1]

        σ, ε, M = 3.4, 120.0, 39.948                # Å, K, g/mol
        eos_si = LJRef(["A"]; userlocations=(; sigma=[σ], epsilon=[ε], Mw=[M], segment=[1.]))
        es_si = ESFramework("A", eos_si; userlocations=coeffs, σε_ce=((σ, ε),))
        σm, εJ, m = σ*1e-10, ε*ES.kB, M*1e-3/ES.NA
        V, TK = ES.NA*σm^3/ϱ, T*ε
        @test ES.VT_viscosity(es_si, V, TK) / (sqrt(m*εJ)/σm^2) ≈ Y[1] rtol=1e-10
        @test ES.VT_thermal_conductivity(es_si, V, TK) / (ES.kB/σm^2*sqrt(εJ/m)) ≈ Y[2] rtol=1e-10
        @test ES.VT_self_diffusion_coefficient(es_si, V, TK)[1] / (σm*sqrt(εJ/m)) ≈ Y[3] rtol=1e-10
    end

    @testset "Fit" begin
        es = ESFramework("LJ", eos; userlocations=coeffs)
        T = [1.5, 2.0, 3.0, 4.0]
        ϱ = [0.3, 0.5, 0.7, 0.8]
        η = ES.VT_viscosity.(es, 1 ./ ϱ, T)
        fit = ESFramework("LJ", eos, [ViscosityData(T, nothing, ϱ, η, :unknown)])
        @test fit[Viscosity()].α1[1] ≈ 1.5 rtol=1e-6
        @test fit[Viscosity()].α2[1] ≈ 0.3 rtol=1e-6
    end
end
