@testset "Framework" begin
    @testset "Fit Pure" begin
        # Load eos and data
        eos = PCSAFT("n-butane")
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
        model = FrameworkModel(eos, data_pure; opts=fit_opts)
        (α_η, α_λ, α_D) = [model[prop_i].α for prop_i in [Viscosity(), ThermalConductivity(), SelfDiffusionCoefficient()]]

        @test all(isapprox(α_η, [0.0; -9.490597; 9.747817; -1.279090; 0.3666153;;]; atol=1e-5))
        @test all(isapprox(α_λ, [3.532197; 116.820461; -99.444029; 24.245996; 0.587130;;]; atol=1e-5))
        @test all(isapprox(α_D, [0.0; 0.0; 0.0; -3.573438; -0.99224;;]; atol=1e-5))
    end
    
    @testset "Fit Inf. Dilution" begin
        # Load eos and data
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
end