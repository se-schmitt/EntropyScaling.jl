
@testset "GCES" begin

    @testset "Pure Viscosity" begin
        # Comparison to table S5 of https://doi.org/10.1007/s10765-022-03096-9
        # acetone
        

        #alcohols

        eos_model = HomogcPCPSAFT(list)
        model = GCESModel(eos_model, list)

        @test viscosity(model_ace,1.3922e4,281.005)/1e-6 ≈ 371.6 rtol=1e-4

        # Test database
        model_ace_db = RefpropRESModel("acetone")
        @test viscosity(model_ace,1e5,300.) ≈ viscosity(model_ace_db,1e5,300.)

        # Martinek et al. (2025)
        model_r11_martinek = RefpropRESModel("R11"; ηref="Martinek et al. (2025)")
        @test viscosity(model_r11_martinek, 9.32e2, 209.08)*1e6 ≈ 1388.935 rtol=1e-4
    end

end



list = [("ethanol", ["CH3" => 1, "CH2" => 1])]#[get_groups_from_name("propanol", gcPCPSAFTGroups)]
println("..")

eos_model = HomogcPCPSAFT(list)
model = GCESModel(eos_model, list)


η = viscosity(model, 1e5, 200.) #in properties
