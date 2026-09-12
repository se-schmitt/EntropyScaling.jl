@testitem "Chapman-Enskog" begin
    @testset "Model" begin
        # Methane 
        σ, ε, Mw = 3.758, 148.6, 16.043             # Poling et al.
        model = ChapmanEnskog("methane"; userlocations=(;sigma=σ, epsilon=ε, Mw))
        @test viscosity(model, 200.)/1e-6 ≈ 7.6848 rtol=1e-2                # NIST value
        model_db = ChapmanEnskog("methane")
        @test viscosity(model, NaN, 300.) ≈ viscosity(model_db, NaN, 300.)

        # Argon
        σ, ε, Mw = 3.350, 143.2, 39.948             # Refprop
        model_Neufeld = ChapmanEnskog("argon"; userlocations=(;sigma=σ, epsilon=ε, Mw), collision_integral=EntropyScaling.Neufeld())
        # @test viscosity(model_Neufeld, 200.)/1e-6 ≈ 7.6848 rtol=1e-2    
        #TODO sth wrong with Neufeld
    end
    
    @testset "Collision Integrals" begin
        @test EntropyScaling.Ω_22(0.9) ≈ 1.68262401248626
        @test EntropyScaling.Ω_11(0.9) ≈ 1.517529517136435
    end
end