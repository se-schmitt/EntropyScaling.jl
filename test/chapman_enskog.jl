@testset "Chapman-Enskog" begin
    @testset "Model" begin
        # Methane 
        σ, ε, Mw = 3.758e-10, 148.6*EntropyScaling.kB, 16.04e-3         # Poling et al.
        model = ChapmanEnskogModel("methane",σ,ε,Mw)
        @test viscosity(model, 200.)/1e-6 ≈ 7.6848 rtol=1e-2            # NIST value

        # Argon
        σ, ε, Mw = 0.3350e-9, 143.2*EntropyScaling.kB, 39.948e-3        # Refprop
        model_Neufeld = ChapmanEnskogModel("argon",σ,ε,Mw;collision_integral=EntropyScaling.Neufeld())
        # @test viscosity(model_Neufeld, 200.)/1e-6 ≈ 7.6848 rtol=1e-2    
        #TODO sth wrong with Neufeld
    end
    
    @testset "Collision Integrals" begin
        @test EntropyScaling.Ω_22(0.9) ≈ 1.68262401248626
        @test EntropyScaling.Ω_11(0.9) ≈ 1.517529517136435
    end
end