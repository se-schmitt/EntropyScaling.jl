@testitem "SymbolicsExt" begin    
    using Clapeyron
    using CoolProp
    using Symbolics

    @variables p, T

    @testset "framework model" begin 
        model_FW = ESFramework("n-butane", PCSAFT("n-butane"))

        @test viscosity(model_FW, p, T) isa Num
        @test self_diffusion_coefficient(model_FW, p, T) isa Num
    end

    @testset "RefpropRES model" begin
        model_RP = RefpropRES("r134a")

        @test viscosity(model_RP, p, T) isa Num
        @test thermal_conductivity(model_RP, p, T) isa Num
    end
    
    @testset "CE model" begin
        model_CE = ChapmanEnskog("decane")

        @test viscosity(model_CE, T) isa Num
    end
end