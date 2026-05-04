@testitem "Framework" begin
    using Clapeyron
    using DelimitedFiles

    @testset "Fit Pure" begin
        # EOS models and data
        eos_model = PCSAFT("n-butane")
        (raw_η, raw_λ, raw_D) = [readdlm("data/exp_$name.csv",',';skipstart=1) for name in ["eta", "lambda", "D"]]

        # Transport properties pure butane
        data_pure = [
            ViscosityData(raw_η[:,1], nothing, raw_η[:,2], raw_η[:,3], :unknown)
            ThermalConductivityData(raw_λ[:,1], nothing, raw_λ[:,2], raw_λ[:,3], :unknown)
            SelfDiffusionCoefficientData(raw_D[:,1], nothing, raw_D[:,2], raw_D[:,3], :unknown)
        ]
        model = ESFramework("n-butane", eos_model, data_pure; tofit=Dict(SelfDiffusionCoefficient()=>[:α2,:α3]))
        _props = [Viscosity(), ThermalConductivity(), SelfDiffusionCoefficient()]
        (α_η, α_λ, α_D) = [first.(getproperty.(model[_p],[:α0,:αln,:α1,:α2,:α3])) for _p in _props]

        α_η_ref = [0.000000e+00;-9.490597e+00; 9.747817e+00;-1.279090e+00; 3.666153e-01;;]
        α_λ_ref = [3.532197e+00; 1.168205e+02;-9.944405e+01; 2.424600e+01; 5.871304e-01;;]
        α_D_ref = [0.000000e+00; 0.000000e+00; 0.000000e+00;-3.573438e+00;-9.922401e-01;;]
        @test isapprox(α_η, α_η_ref; rtol=1e-5)
        @test isapprox(α_λ, α_λ_ref; rtol=1e-5)
        @test isapprox(α_D, α_D_ref; rtol=1e-5)
    end
    
    @testset "Fit Inf. Dilution" begin
        # EOS models and data
        system = PCSAFT(["hexane", "toluene"])
        raw = readdlm("data/exp_D_inf.csv",',';skipstart=1) 

        data_inf = [
            TransportPropertyData(InfDiffusionCoefficient(1 => 2), raw[:,1], nothing, raw[:,2], raw[:,3], :unknown)
        ]
        model = ESFramework(["hexane", "toluene"], system, data_inf; tofit=Dict(InfDiffusionCoefficient(1=>2)=>[:α2,:α3]))

        α_D_inf = first.(getproperty.(model[InfDiffusionCoefficient(1 => 2)],[:α0,:αln,:α1,:α2,:α3]))
        α_D_inf_ref = [0.0; 0.0; 0.0;-2.604944e+00;-1.567851e+00;;]
        @test all(isapprox(α_D_inf, α_D_inf_ref; rtol=1e-5))
    end

    @testset "Call" begin
        @testset "Pure" begin
            # EOS model
            eos_model = PCSAFT("n-butane")
            
            # ES model 
            model = ESFramework("n-butane", eos_model; 
                userlocations = Dict(
                    Viscosity() => (; αln=[-14.165], α1=[13.97], α2=[-2.382], α3=[0.501]),
                    ThermalConductivity() => (; α0=[3.962], αln=[98.222], α1=[-82.974], α2=[20.079], α3=[1.073]),
                    SelfDiffusionCoefficient() => (; α2=[-3.507], α3=[-0.997]),
            ))

            @test viscosity(model, 37.21e6, 323.) ≈ 1.921922e-4 rtol=1e-5
            @test thermal_conductivity(model, 37.21e6, 323.) ≈ 1.199070e-1 rtol=1e-5
            @test self_diffusion_coefficient(model, 37.21e6, 323.) ≈ 6.645588e-9 rtol=1e-5
        end

        @testset "Mixtures" begin
            # ES models
            model_1 = ESFramework(["benzene","hexane"], PCSAFT(["benzene","hexane"]);
                userlocations=Dict(
                    ThermalConductivity() => (; α0=[6.492054,9.756965], α2=[1.985553,1.657211], α3=[3.126438,5.058428])
            ))
            model_2 = ESFramework(["hexane", "dodecane"], PCSAFT(["hexane", "dodecane"]);
                userlocations=Dict(
                    DiffusionCoefficient() => (; α2=[-2.5414 -4.0610;-2.5463 -4.6610], α3=[-1.9186 -3.4267;-1.8070 -3.2021]),
            ))

            @test thermal_conductivity(model_1, 0.1e6, 294.7, [.5,.5]) ≈ 1.338512e-1 rtol=1e-5
            @test isapprox(self_diffusion_coefficient(model_2, 0.1e6, 298.15, [.25,.75]),[1.746913e-9,1.097029e-9]; rtol=1e-5)
            @test MS_diffusion_coefficient(model_2, 0.1e6, 298.15, [.6,.4])[1,2] ≈ 1.912959e-9 rtol=1e-5
            @test self_diffusion_coefficient(model_2, 0.1e6, 298.15, [0.,1.])[1] ≈ MS_diffusion_coefficient(model_2, 0.1e6, 298.15, [0.,1.])[1,2]
            @test fick_diffusion_coefficient(model_2, 0.1e6, 298.15, [.6,.4])[1,1] ≈ 1.967921e-9 rtol=1e-5
        end
    end
end