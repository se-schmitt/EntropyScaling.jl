
@testset "RefpropRES" begin
    @testset "Pure" begin
        # Comparison to table S5 of https://doi.org/10.1007/s10765-022-03096-9
        # acetone
        ni_ace = [-0.501124;1.194414;-0.521023;0.079149;;]
        ξ_ace = [1.]
        σ_ace = [0.467e-9]
        ε_ace = [519*EntropyScaling.kB]

        model_ace = RefpropRESModel(SingleFluid("acetone"),Dict(
            Viscosity() => [ni_ace, ξ_ace, σ_ace, ε_ace]
        ))
        @test viscosity(model_ace,1.3922e4,281.005)/1e-6 ≈ 371.6 rtol=1e-4

        # cyclopentane
        ni_cyc = [-0.448046;1.012681;-0.381869;0.054674;;]
        ξ_cyc = [0.9492]
        σ_cyc = [0.518e-9]
        ε_cyc = [406.33*EntropyScaling.kB]

        model_cyc = RefpropRESModel(SingleFluid("cyclopentane"),Dict(
            Viscosity() => [ni_cyc, ξ_cyc, σ_cyc, ε_cyc]
        ))
        @test viscosity(model_cyc,1.4244e4,273.15)/1e-6 ≈ 567.0 rtol=1e-4

        # Test database
        model_ace_db = RefpropRESModel("acetone")
        @test viscosity(model_ace,1e5,300.) ≈ viscosity(model_ace_db,1e5,300.)

        model_cyc_db = RefpropRESModel("cyclopentane")
        @test viscosity(model_cyc,1e5,300.) ≈ viscosity(model_cyc_db,1e5,300.)
    end

    @testset "Mixtures" begin
        # Comparison to table S7 of https://doi.org/10.1007/s10765-022-03096-9
        # decane + methane
        ni_dec = [-0.435268;0.876185;-0.314288;0.039133;;]
        ni_met = [0.140367;-0.046685;0.245800;-0.050594;;]
        ni_mix = [ni_dec ni_met]
        ξ_mix = [1., 1.]
        σ_mix = [0.686e-9, 0.367e-9]
        ε_mix = [490.51, 174].*EntropyScaling.kB

        model_mix = RefpropRESModel(MultiFluid(["decane","methane"]),Dict(
            Viscosity() => [ni_mix, ξ_mix, σ_mix, ε_mix]
        ))
        @test viscosity(model_mix,3.565407e6,289.996,[0.8680,1-0.8680])/1e-6 ≈ 925.02 rtol = 1e-5

        # Test database
        model_mix_db = RefpropRESModel(["decane","methane"])
        @test viscosity(model_mix, 1e5, 300., [.5,.5]) ≈ viscosity(model_mix_db, 1e5, 300., [.5,.5])
    end 
end