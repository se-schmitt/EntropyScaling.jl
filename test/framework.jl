@testset "Framework" begin
    @testset "Fit Pure" begin
        # EOS models and data
        pure = PCSAFT("n-butane")
        (raw_η, raw_λ, raw_D) = [readdlm("data/exp_$name.csv",',';skipstart=1) for name in ["eta", "lambda", "D"]]

        # Transport properties pure butane
        data_pure = [
            ViscosityData(raw_η[:,1], [], raw_η[:,2]./0.05812, raw_η[:,3], :unknown)
            ThermalConductivityData(raw_λ[:,1], [], raw_λ[:,2]./0.05812, raw_λ[:,3], :unknown)
            SelfDiffusionCoefficientData(raw_D[:,1], [], raw_D[:,2]./0.05812, raw_D[:,3], :unknown)
        ]
        fit_opts = FitOptions(;
            what_fit=Dict(
                ThermalConductivity()=>ones(Bool,5), 
                SelfDiffusionCoefficient()=>Bool[0,0,0,1,1])
        )
        model = FrameworkModel(pure, data_pure; opts=fit_opts)
        (α_η, α_λ, α_D) = [model[prop_i].α for prop_i in [Viscosity(), ThermalConductivity(), SelfDiffusionCoefficient()]]

        @test all(isapprox(α_η, [0.0; -9.490597; 9.747817; -1.279090; 0.3666153;;]; atol=1e-5))
        @test all(isapprox(α_λ, [3.532197; 116.820461; -99.444029; 24.245996; 0.587130;;]; atol=1e-5))
        @test all(isapprox(α_D, [0.0; 0.0; 0.0; -3.573438; -0.99224;;]; atol=1e-5))
    end
    
    @testset "Fit Inf. Dilution" begin
        # EOS models and data
        solvent = PCSAFT("toluene")
        solute = PCSAFT("hexane")
        raw = readdlm("data/exp_D_inf.csv",',';skipstart=1) 

        data_inf = [
            InfDiffusionCoefficientData(raw[:,1], [], raw[:,2]./0.09214, raw[:,3], :unknown)
        ]
        fit_opts = FitOptions(;what_fit=Dict(InfDiffusionCoefficient()=>Bool[0,0,0,1,1]))
        model = FrameworkModel(solvent, data_inf; opts=fit_opts, solute=solute)

        @test all(isapprox(model[InfDiffusionCoefficient()].α, [0.0; 0.0; 0.0; -2.604948; -1.567851;;]; atol=1e-5))
    end

    @testset "Call" begin
        @testset "Pure" begin
            # EOS model
            pure = PCSAFT("n-butane")
            
            # ES model 
            model = FrameworkModel(pure,Dict(
                Viscosity() => [0.;-14.165;13.97;-2.382;0.501;;],
                ThermalConductivity() => [3.962;98.222;-82.974;20.079;1.073;;],
                SelfDiffusionCoefficient() => [0.;0.;0.;-3.507;-0.997;;]
            ))

            @test viscosity(model, pure, 37.21e6, 323.) ≈ 1.921922e-4 rtol=1e-4
            @test thermal_conductivity(model, pure, 37.21e6, 323.) ≈ 1.199070e-1 rtol=1e-4
            @test self_diffusion_coefficient(model, pure, 37.21e6, 323.) ≈ 6.645588e-9 rtol=1e-4
        end

        @testset "Mixtures" begin
            # EOS models     
            mix1 = PCSAFT(["benzene","hexane"])
            user_mix2 = (; dipole=[2.88,0.0], epsilon=[232.99,271.63], sigma=[3.2742,3.4709], Mw=[58.08,119.38], segment=[2.7447,2.5038])
            mix2 = PCPSAFT(["acetone", "chloroform"]; userlocations=user_mix2)
            mix3 = PCSAFT(["toluene","heptane"])

            # ES models
            model_1 = FrameworkModel(mix1, Dict(
                ThermalConductivity() => hcat(
                    [6.492054, 0.0, 0.0, 1.985553, 3.126438],
                    [9.756965, 0.0, 0.0, 1.657211, 5.058428]
                )
            ))

            @test thermal_conductivity(model_1, mix1, )
        end
        # α_λ_mix1 = [
        #     [6.492054425320112,0.0,0.0,1.9855528414772259,3.12643833453098],
        #     [9.75696539716926,0.0,0.0,1.6572105176108498,5.058427863917941]
        # ]
        # α_D_mix2 = [
        #     [0.,0.,0.,1.1301081958e+01,-7.7638664176e+00],
        #     [0.,0.,0.,-3.0212807487e+00,-1.8500112748e+00]
        # ]
        # α_Ð_mix3 = [
        #     [0.,0.,0.,1.1606624185e+01,-7.4861787912e+00],
        #     [0.,0.,0.,1.8829820233e+01,-1.3186835311e+01],
        # ]
    end
end