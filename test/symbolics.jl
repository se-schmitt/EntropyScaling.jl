@testset "SymbolicsExt" begin    
    @variables p, T

    @testset "framework model" begin 
        model_FW = FrameworkModel(PCSAFT("n-butane"),Dict(
            DynamicViscosity() => [[0.;-14.165;13.97;-2.382;0.501;;]],
            SelfDiffusionCoefficient() => [[0.;0.;0.;-3.507;-0.997;;]]
        ))

        @test dynamic_viscosity(model_FW, p, T) isa Num
        @test self_diffusion_coefficient(model_FW, p, T) isa Num
    end

    @testset "RefpropRES model" begin
        model_RP = RefpropRESModel("r134a")

        @test dynamic_viscosity(model_RP, p, T) isa Num
        @test thermal_conductivity(model_RP, p, T) isa Num
    end
    
    @testset "CE model" begin
        model_CE = ChapmanEnskogModel("decane")

        @test dynamic_viscosity(model_CE, T) isa Num
    end
end