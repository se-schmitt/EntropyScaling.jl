@testset "Chapman-Enskog" begin
    @testset "LJ Collision Int" begin
        @test EntropyScaling.Ω_22(0.9) ≈ 1.68262401248626
        @test EntropyScaling.Ω_11(0.9) ≈ 1.517529517136435
    end
end