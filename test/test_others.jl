@testitem "Others" begin
    @testset "PolynomialDiluteGas" begin
        model_but = PolynomialDiluteGas("butane")
        model = PolynomialDiluteGas(["butane", "argon"])

        T = 300.0
        x = [0.25, 0.75]

        @test viscosity(model, NaN, T, [1.0, 0.0]) ≈ viscosity(model_but, NaN, T)
        @test viscosity(model, NaN, T, x) ≈ 1.5663924445319993e-5 rtol=1e-8
    end
end