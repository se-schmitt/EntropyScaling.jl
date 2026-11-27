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
    
    @testset "Maxwell-Stefan diffusion" begin
        σ, ε, M = [3.758,3.542,3.798], [148.6,93.3,71.4], [16.043,39.948,28.014]
        model = ChapmanEnskog(["methane","argon","nitrogen"]; userlocations=(; sigma=σ, epsilon=ε, Mw=M))
        T, p = 300.0, 101325.0
        ϱD = MS_diffusion_coefficient(model, p, T, [.3,.3,.4])
        for (i,j) in ((1,2),(1,3),(2,3))
            σij, εij, Mij = (σ[i]+σ[j])/2, sqrt(ε[i]*ε[j]), 2/(1/M[i]+1/M[j])
            D_hirschfelder = 0.0018583*sqrt(T^3*(1/M[i]+1/M[j]))/(σij^2*EntropyScaling.Ω_11(T/εij))*1e-4
            @test ϱD[i,j] / (p/(EntropyScaling.kB*T) * Mij*1e-3/EntropyScaling.NA) ≈ D_hirschfelder rtol=1e-3
        end
    end

    @testset "Collision Integrals" begin
        @test EntropyScaling.Ω_22(0.9) ≈ 1.68262401248626
        @test EntropyScaling.Ω_11(0.9) ≈ 1.517529517136435
    end
end