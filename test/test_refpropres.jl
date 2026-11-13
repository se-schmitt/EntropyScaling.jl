
@testitem "RefpropRES" begin
    using Clapeyron
    using CoolProp
    using Downloads

    const url_refprop = "https://raw.githubusercontent.com/usnistgov/fastchebpure/50af5c154a113ac27a2c0a1c3538bc4f43a73a66/teqp_REFPROP10/dev/fluids/"

    @testset "Pure Viscosity" begin
        # Comparison to table S5 of https://doi.org/10.1007/s10765-022-03096-9
        # acetone
        model_ace = RefpropRES("acetone", SingleFluid("acetone"); 
            userlocations=Dict(
                Viscosity() => (;n1=[-0.501124], n2=[1.194414], n3=[-0.521023], n4=[0.079149], ξ=[1.])
            ),
            ce_userlocations=(;sigma=[4.67], epsilon=[519.], Mw=[58.08])
        )
        @test viscosity(model_ace,1.3922e4,281.005)/1e-6 ≈ 371.6 rtol=1e-4

        # cyclopentane
        model_cyc = RefpropRES("cyclopentane", SingleFluid("cyclopentane");
            userlocations=Dict(
                Viscosity() => (;n1=[-0.448046], n2=[1.012681], n3=[-0.381869], n4=[0.054674], ξ=[0.9492])
            ),
            ce_userlocations=(;sigma=[5.18], epsilon=[406.33], Mw=[70.1329])
        )
        @test viscosity(model_cyc,1.4244e4,273.15)/1e-6 ≈ 567.0 rtol=1e-4

        # Test database
        model_ace_db = RefpropRES("acetone")
        @test viscosity(model_ace,1e5,300.) ≈ viscosity(model_ace_db,1e5,300.)

        model_cyc_db = RefpropRES("cyclopentane")
        @test viscosity(model_cyc,1e5,300.) ≈ viscosity(model_cyc_db,1e5,300.)

        buta13_path = joinpath(url_refprop,"13BUTADIENE.json")
        buta13_json = Downloads.download(buta13_path) |> read |> String
        model_13but = RefpropRES("1,3-butadiene", SingleFluid(buta13_json; coolprop_userlocations=false))
        @test viscosity(model_13but,9.259e3,220.)/1e-6 ≈ 424.7 rtol=1e-4

        # Martinek et al. (2025)
        model_r11_martinek = RefpropRES2025("R11")
        @test viscosity(model_r11_martinek, 9.32e2, 209.08)*1e6 ≈ 1388.935 rtol=1e-4
    end

    @testset "Mixtures Viscosity" begin
        # Comparison to table S7 of https://doi.org/10.1007/s10765-022-03096-9
        # decane + methane
        model_mix = RefpropRES(["decane","methane"], MultiFluid(["decane","methane"]);
            userlocations=Dict(
                Viscosity() => (; n1=[-0.435268,0.140367], n2=[0.876185,-0.046685], n3=[-0.314288,0.245800], n4=[0.039133,-0.050594], ξ=[1.,1.]),
            ),
            ce_userlocations=(;sigma=[6.86,3.67], epsilon=[490.51, 174], Mw=[142.282,16.0428])
        )
        @test viscosity(model_mix,3.565407e6,289.996,[0.8680,1-0.8680])/1e-6 ≈ 925.02 rtol = 1e-5

        # Test database
        model_mix_db = RefpropRES(["decane","methane"])
        @test viscosity(model_mix, 1e5, 300., [.5,.5]) ≈ viscosity(model_mix_db, 1e5, 300., [.5,.5])
    end 

    @testset "Pure Th. Cond." begin
        model_hex = RefpropRES("n-hexane")
        model_r13 = RefpropRES("R13")

        # Test against Yang et al. (2021) (10.1021/acs.iecr.1c02154)
        @test thermal_conductivity(model_hex, 0.101e6, 183.16) ≈ 0.1571 atol=1e-4
        @test thermal_conductivity(model_r13, 0.109e6, 193.10) ≈ 0.0886 atol=1e-4

        # Test against Li et al. (2024) (10.1021/acs.iecr.4c02946)
        model_amm = RefpropRES2025("ammonia")
        model_D2 = RefpropRES2025("deuterium")
        model_R1233ZDE = RefpropRES2025("R1233ZDE")

        @test thermal_conductivity(model_amm, 23050e3, 206.195) ≈ 0.6495 atol=1e-5
        @test thermal_conductivity(model_D2, 101.1e3, 338.15) ≈ 0.14403 atol=1e-4
        @test thermal_conductivity(model_R1233ZDE, 7100e3, 242.54) ≈ 0.10332 atol=1e-5
    end

    @testset "Mix Th. Cond." begin
        model_mix = RefpropRES(["R134a","R1234yf"])

        # Test against Yang et al. (2021) (10.1021/acs.iecr.1c02154)
        @test thermal_conductivity(model_mix, 1.01e6, 254.86, [0.504,0.496]) ≈ 0.0879 atol=1e-3

        # Test against Li et al. (2024) (10.1021/acs.iecr.4c02946)
        model_mix = RefpropRES2025(["ethanol","water"])

        @test thermal_conductivity(model_mix, 101e3, 284.143, [0.28111,0.71889]) ≈ 0.28867 atol=1e-4
    end
end