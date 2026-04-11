@testitem "misc" begin
    @testset "ConstantModel" begin
        vec = [
            viscosity => Viscosity(), 
            thermal_conductivity => ThermalConductivity(), 
            self_diffusion_coefficient => SelfDiffusionCoefficient(),
            MS_diffusion_coefficient => MaxwellStefanDiffusionCoefficient(), 
            fick_diffusion_coefficient => FickDiffusionCoefficient(),
            inf_diffusion_coefficient => InfDiffusionCoefficient(),
        ]
        for (fun,prop) in vec
            v = rand()
            model = ConstantModel(prop, v)
            @test fun(model, nothing, nothing) == v
        end
    end
end