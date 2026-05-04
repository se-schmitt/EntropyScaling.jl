@testitem "GCES" begin
    using Clapeyron

    @testset "pure vis" begin
        # ethane
        model = GCES(("ethane", ["CH3" => 2]))
        T, p = 250.0, 10e6
        @test viscosity(model, p, T) ≈ 0.087084e-3 rtol=1e-3

        # hexane
        comp  = ("hexane", ["CH3" => 2, "CH2" => 4])
        eos   = HomogcPCPSAFT(comp)
        model = GCES(comp, eos)
        T, p = 450.0, 4e6
        @test viscosity(model, p, T) ≈ 0.088013e-3 rtol=1e-3

        # benzene
        model = GCES(("benzene", ["aCH" => 6]))
        T, p = 400.0, 10e6
        @test viscosity(model, p, T) ≈ 0.234569e-3 rtol=2e-3

        # cyclopentane
        comp  = ("cyclopentane", ["cCH2_pent" => 5])
        eos   = HomogcPCPSAFT([comp])
        model = GCES(comp, eos)
        T, p = 350.0, 10e6
        @test viscosity(model, p, T) ≈ 0.286097e-3 rtol=2e-3

        # cyclohexane
        model = GCES("cyclohexane")
        T, p = 400.0, 10e6
        @test viscosity(model, p, T) ≈ 0.300706e-3 rtol=2e-3
    end

    @testset "mix vis" begin
        model = GCES(["hexane", "benzene"])
        T, p, x = 300.0, 1e6, [.25,.75]
        @test viscosity(model, p, T, x) ≈ 0.0004404487489691302 rtol=1e-8
    end
end