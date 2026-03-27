@testset "GCES" begin
    # ethane
    comp  = ("ethane", ["CH3" => 2])
    eos   = HomogcPCPSAFT([comp])
    model = GCESModel(eos, [comp])
    T = 250.0; p = 10e6
    @test Clapeyron.mass_density(eos, p, T) ≈ 555.563524 rtol=1e-4
    @test viscosity(model, p, T) ≈ 0.087084e-3 rtol=1e-3

    # hexane
    comp  = ("hexane", ["CH3" => 2, "CH2" => 4])
    eos   = HomogcPCPSAFT([comp])
    model = GCESModel(eos, [comp])
    T = 450.0; p = 4e6
    @test Clapeyron.mass_density(eos, p, T) ≈ 484.770693 rtol=1e-4
    @test viscosity(model, p, T) ≈ 0.088013e-3 rtol=1e-3

    # benzene
    comp  = ("benzene", ["aCH" => 6])
    eos   = HomogcPCPSAFT([comp])
    model = GCESModel(eos, [comp])
    T = 400.0; p = 10e6
    @test Clapeyron.mass_density(eos, p, T) ≈ 691.279465 rtol=1e-4
    @test viscosity(model, p, T) ≈ 0.234569e-3 rtol=2e-3

    # cyclopentane
    comp  = ("cyclopentane", ["cCH2_pent" => 5])
    eos   = HomogcPCPSAFT([comp])
    model = GCESModel(eos, [comp])
    T = 350.0; p = 10e6
    @test Clapeyron.mass_density(eos, p, T) ≈ 701.666163 rtol=1e-4
    @test viscosity(model, p, T) ≈ 0.286097e-3 rtol=2e-3

    # cyclohexane
    comp  = ("cyclohexane", ["cCH2_hex" => 6])
    eos   = HomogcPCPSAFT([comp])
    model = GCESModel(eos, [comp])
    T = 400.0; p = 10e6
    @test Clapeyron.mass_density(eos, p, T) ≈ 701.151517 rtol=1e-4
    @test viscosity(model, p, T) ≈ 0.300706e-3 rtol=2e-3
end