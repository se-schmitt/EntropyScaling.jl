@testset "Framework" begin
    @testset "Fit Pure" begin
        # EOS models and data
        eos_model = PCSAFT("n-butane")
        (raw_η, raw_λ, raw_D) = [readdlm("data/exp_$name.csv",',';skipstart=1) for name in ["eta", "lambda", "D"]]

        # Transport properties pure butane
        data_pure = [
            ViscosityData(raw_η[:,1], [], raw_η[:,2], raw_η[:,3], :unknown)
            ThermalConductivityData(raw_λ[:,1], [], raw_λ[:,2], raw_λ[:,3], :unknown)
            SelfDiffusionCoefficientData(raw_D[:,1], [], raw_D[:,2]
            , raw_D[:,3], :unknown)
        ]
        fit_opts = FitOptions(;
            what_fit=Dict(
                ThermalConductivity()=>ones(Bool,5), 
                SelfDiffusionCoefficient()=>Bool[0,0,0,1,1])
        )
        model = FrameworkModel(eos_model, data_pure; opts=fit_opts)
        (α_η, α_λ, α_D) = [model[prop_i].α for prop_i in [Viscosity(), ThermalConductivity(), SelfDiffusionCoefficient()]]

        α_η_ref = [0.000000e+00;-9.490597e+00; 9.747817e+00;-1.279090e+00; 3.666153e-01;;]
        α_λ_ref = [3.532197e+00; 1.168205e+02;-9.944405e+01; 2.424600e+01; 5.871304e-01;;]
        α_D_ref = [0.000000e+00; 0.000000e+00; 0.000000e+00;-3.573438e+00;-9.922401e-01;;]
        @test isapprox(α_η, α_η_ref; rtol=1e-5)
        @test isapprox(α_λ, α_λ_ref; rtol=1e-5)
        @test isapprox(α_D, α_D_ref; rtol=1e-5)
    end
    
    @testset "Fit Inf. Dilution" begin
        # EOS models and data
        solvent = PCSAFT("toluene")
        solute = PCSAFT("hexane")
        raw = readdlm("data/exp_D_inf.csv",',';skipstart=1) 

        data_inf = [
            InfDiffusionCoefficientData(raw[:,1], [], raw[:,2], raw[:,3], :unknown)
        ]
        fit_opts = FitOptions(;what_fit=Dict(InfDiffusionCoefficient()=>Bool[0,0,0,1,1]))
        model = FrameworkModel(solvent, data_inf; opts=fit_opts, solute=solute)

        α_D_inf_ref = [0.0; 0.0; 0.0;-2.604944e+00;-1.567851e+00;;]
        @test all(isapprox(model[InfDiffusionCoefficient()].α, α_D_inf_ref; rtol=1e-5))
    end

    @testset "Call" begin
        @testset "Pure" begin
            # EOS model
            eos_model = PCSAFT("n-butane")
            
            # ES model 
            model = FrameworkModel(eos_model,Dict(
                Viscosity() => [0.;-14.165;13.97;-2.382;0.501;;],
                ThermalConductivity() => [3.962;98.222;-82.974;20.079;1.073;;],
                SelfDiffusionCoefficient() => [0.;0.;0.;-3.507;-0.997;;]
            ))

            @test viscosity(model, 37.21e6, 323.) ≈ 1.921922e-4 rtol=1e-5
            @test thermal_conductivity(model, 37.21e6, 323.) ≈ 1.199070e-1 rtol=1e-5
            @test self_diffusion_coefficient(model, 37.21e6, 323.) ≈ 6.645588e-9 rtol=1e-5
        end

        @testset "Mixtures" begin
            # ES models
            model_1 = FrameworkModel(PCSAFT(["benzene","hexane"]), Dict(
                ThermalConductivity() => hcat(
                    [6.492054, 0.0, 0.0, 1.985553, 3.126438],
                    [9.756965, 0.0, 0.0, 1.657211, 5.058428]),
            ))
            model_2 = FrameworkModel(PCSAFT(["hexane", "dodecane"]), Dict(
                SelfDiffusionCoefficient() => hcat(
                    [0.0, 0.0, 0.0, -2.5414, -1.9186],
                    [0.0, 0.0, 0.0, -4.6610, -3.2021]
                ),
                InfDiffusionCoefficient() => hcat(
                    [0.0, 0.0, 0.0, -2.5463, -1.8070],
                    [0.0, 0.0, 0.0, -4.0610, -3.4267]),
            ))

            @test thermal_conductivity(model_1, 0.1e6, 294.7, [.5,.5]) ≈ 1.338512e-1 rtol=1e-5
            @test isapprox(self_diffusion_coefficient(model_2, 0.1e6, 298.15, [.25,.75]),[1.746913e-9,1.097029e-9]; rtol=1e-5)
            @test MS_diffusion_coefficient(model_2, 0.1e6, 298.15, [.6,.4]) ≈ 1.912959e-9 rtol=1e-5
            @test self_diffusion_coefficient(model_2, 0.1e6, 298.15, [0.,1.])[1] ≈ MS_diffusion_coefficient(model_2, 0.1e6, 298.15, [0.,1.])
        end
    end
end