@testset "UnitfulExt" begin

    @testset "Data Constructors" begin
        ηdatu = ViscosityData([300.0u"K"], [1u"bar"], [1u"mPa*s"])
        ηdatU = TransportPropertyData([300.0u"K"], [1u"bar"], [1u"cP"])
        ηdat = ViscosityData([300.], [1e5], nothing, [1e-3])
        @test all([getfield(ηdatu,sym) == getfield(ηdatU,sym) == getfield(ηdat,sym) for sym in [:T,:p,:Y]])

        λdatu = ThermalConductivityData([25u"°C"], [10u"mol/l"], [0.15u"W/(m*K)"])
        λdatU = TransportPropertyData([25u"°C"], [10u"mol/l"], [0.15u"W/(m*K)"])
        λdat = ThermalConductivityData([298.15], nothing, [10000.], [0.15])
        @test all([getfield(λdatu,sym) == getfield(λdatU,sym) == getfield(λdat,sym) for sym in [:T,:ϱ,:Y]])

        Ddatu = SelfDiffusionCoefficientData([80.33u"°F"], [1u"bar"], [1e-5u"cm^2/s"])
        Ddat = SelfDiffusionCoefficientData([300.], [1e5], nothing, [1e-9])
        @test all([getfield(Ddatu,sym) == getfield(Ddat,sym) for sym in [:T,:p,:Y]])

        Dinfdatu = InfDiffusionCoefficientData([80.33u"°F"], [10u"mol/dm^3"], [1e-5u"cm^2/s"])
        Dinfdat = InfDiffusionCoefficientData([300.], nothing, [10000.], [1e-9])
        @test all([getfield(Dinfdatu,sym) == getfield(Dinfdat,sym) for sym in [:T,:ϱ,:Y]])
    end

    @testset "Property functions" begin

        model = FrameworkModel(PCSAFT("n-butane"),Dict(
            Viscosity() => [[0.;-14.165;13.97;-2.382;0.501;;]],
            ThermalConductivity() => [[3.962;98.222;-82.974;20.079;1.073;;]],
            SelfDiffusionCoefficient() => [[0.;0.;0.;-3.507;-0.997;;]]
        ))

        @test viscosity(model, 37.21u"MPa", 323u"K"; output=u"cP").val ≈ 1.921922e-1 rtol=1e-5
        @test thermal_conductivity(model, 372.1u"bar", 49.85u"°C").val ≈ 1.199070e-1 rtol=1e-5
        @test self_diffusion_coefficient(model, 37.21u"MPa", 323u"K"; output=u"cm^2/s").val ≈ 6.645588e-5 rtol=1e-5

        ϱmass = 602.436942723435u"kg/m^3"       #mass_density(model.eos, 37.21u"MPa", 323u"K")
        ϱmol = 10365.398188634463u"mol/m^3"     #molar_density(model.eos, 37.21u"MPa", 323u"K")
        @test viscosity(model, ϱmass, 323u"K").val ≈ 1.921922e-4 rtol=1e-5
        @test thermal_conductivity(model, ϱmol, 49.85u"°C").val ≈ 1.199070e-1 rtol=1e-5
    end
end
