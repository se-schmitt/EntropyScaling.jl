@testitem "Lennard-Jones units" begin
    using Clapeyron
    ES = EntropyScaling
    σ, ε, M = 3.4, 120.0, 39.948                    # LJ unit system [Å, K, g/mol]
    σm, εJ, m = σ*1e-10, ε*ES.kB, M*1e-3/ES.NA
    ηu, λu, Du = sqrt(m*εJ)/σm^2, ES.kB/σm^2*sqrt(εJ/m), σm*sqrt(εJ/m)
    _VT(v, T) = (ES.NA*σm^3*v, T*ε)

    @testset "Chapman-Enskog" begin
        ce = ChapmanEnskog("A"; userlocations=(; sigma=1., epsilon=1., Mw=1.), units=ES.LJ_UNITS)
        for T in (0.8, 2.0, 5.0)
            η = viscosity(ce, NaN, T)
            @test η * ES.Ω(Viscosity(), ce, T; i=1) ≈ 5/16 * sqrt(T/π)
            @test thermal_conductivity(ce, NaN, T) / η ≈ 15/4
        end
    end

    @testset "Pure" begin
        c = Dict(
            Viscosity() => (; α0=0., α1=1.5, α2=0.3, α3=0., αln=0.),
            ThermalConductivity() => (; α0=1., α1=1.2, α2=0.2, α3=0., αln=0.),
            SelfDiffusionCoefficient() => (; α0=zeros(1,1), α1=fill(-1.,1,1), α2=fill(-.2,1,1), α3=zeros(1,1), αln=zeros(1,1)))
        eos = KolafaNezbeda("A"; molecular_units=true, userlocations=(; sigma=[1.], epsilon=[1.], Mw=[1.]))
        eos_si = KolafaNezbeda("A"; userlocations=(; sigma=[σ], epsilon=[ε], Mw=[M]))
        es, es_si = ESFramework("A", eos; userlocations=c), ESFramework("A", eos_si; userlocations=c)
        @test es[Viscosity()].ce.units === ES.LJ_UNITS

        v, T = 1/0.6, 2.0
        V, TK = _VT(v, T)
        η = ES.VT_viscosity(es, v, T)
        @test ES.VT_viscosity(es_si, V, TK)/ηu ≈ η rtol=1e-10
        @test ES.VT_thermal_conductivity(es_si, V, TK)/λu ≈ ES.VT_thermal_conductivity(es, v, T) rtol=1e-10
        @test ES.VT_self_diffusion_coefficient(es_si, V, TK)[1]/Du ≈ ES.VT_self_diffusion_coefficient(es, v, T)[1] rtol=1e-10
        @test viscosity(es, pressure(eos, v, T), T; phase=:liquid) ≈ η

        # fit recovers the parameters
        Ts, ϱs = [1.5, 2.0, 3.0, 4.0], [0.3, 0.5, 0.7, 0.8]
        ηs = ES.VT_viscosity.(es, 1 ./ ϱs, Ts)
        fit = ESFramework("A", eos, [ViscosityData(Ts, nothing, ϱs, ηs, :unknown)])
        @test fit[Viscosity()].α1[1] ≈ 1.5 rtol=1e-6
        @test fit[Viscosity()].α2[1] ≈ 0.3 rtol=1e-6
    end

    @testset "Mixture" begin
        σs, εs, ms = [1., 1.1], [1., .85], [1., 2.]
        c = Dict(
            Viscosity() => (; α0=[0.,0.], α1=[1.5,1.4], α2=[.3,.3], α3=[0.,0.], αln=[0.,0.]),
            ThermalConductivity() => (; α0=[1.,1.], α1=[1.2,1.1], α2=[.2,.2], α3=[0.,0.], αln=[0.,0.]),
            DiffusionCoefficient() => (; α0=zeros(2,2), α1=[-1. -1.1; -.9 -1.05], α2=fill(-.2,2,2), α3=zeros(2,2), αln=zeros(2,2)))
        eos = KolafaNezbeda(["A","B"]; molecular_units=true, userlocations=(; sigma=σs, epsilon=εs, Mw=ms))
        eos_si = KolafaNezbeda(["A","B"]; userlocations=(; sigma=σs.*σ, epsilon=εs.*ε, Mw=ms.*M))
        es, es_si = ESFramework(["A","B"], eos; userlocations=c), ESFramework(["A","B"], eos_si; userlocations=c)

        v, T, z = 1/0.6, 2.0, [.3, .7]
        V, TK = _VT(v, T)
        η = ES.VT_viscosity(es, v, T, z)
        @test ES.VT_viscosity(es_si, V, TK, z)/ηu ≈ η rtol=1e-10
        @test ES.VT_thermal_conductivity(es_si, V, TK, z)/λu ≈ ES.VT_thermal_conductivity(es, v, T, z) rtol=1e-10
        @test ES.VT_self_diffusion_coefficient(es_si, V, TK, z)./Du ≈ ES.VT_self_diffusion_coefficient(es, v, T, z) rtol=1e-10
        @test ES.VT_MS_diffusion_coefficient(es_si, V, TK, z)[1,2]/Du ≈ ES.VT_MS_diffusion_coefficient(es, v, T, z)[1,2] rtol=1e-10
        @test viscosity(es, pressure(eos, v, T, z), T, z; phase=:liquid) ≈ η
    end
end
