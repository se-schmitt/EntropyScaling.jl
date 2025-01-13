@testset "Unitful Extension" begin

    import Unitful: K, °C, bar, atm, Pa, s, mol, m, percent, W
    import Unitful: btu, hr, ft, Ra, cP, μPa, mW, psi, Ra

    Ts = collect(LinRange(10, 80, 15) * °C)
    Ps = [1bar, 2bar, 5bar, 50bar, 1000.0bar]

    fluid = "Toluene"

    ϱ = typeof(1.0 * mol / m^3)[]
    λ = typeof(1.0 * W / m / K)[]
    η = typeof(1.0 * Pa * s)[]

    PT = []
    P_data = typeof(1.0Pa)[]
    T_data = typeof(1.0K)[]
    for T in Ts
        for p in Ps
            push!(ϱ, PropsSI("DMOLAR", "P", p, "T", T, fluid))
            push!(λ, PropsSI("L", "P", p, "T", T, fluid))
            push!(η, PropsSI("VISCOSITY", "P", p, "T", T, fluid))
            push!(PT, (p, T))
            push!(P_data, p)
            push!(T_data, T)
        end
    end


    data_eta_unit_density = ViscosityData(T_data, ϱ, η)
    data_lamb_unit_density = ThermalConductivityData(T_data, ϱ, λ)
    data_eta_unit_pressure = ViscosityData(T_data, P_data, η)
    data_lamb_unit_pressure = ThermalConductivityData(T_data, P_data, λ)

    @testset "Viscosity data generation" begin


        data_eta_ustrip_pressure = ViscosityData(
            ustrip.(T_data .|> K),
            ustrip.(P_data .|> Pa),
            [],
            ustrip.(η),
            :unknown,
        )
        data_eta_ustrip_density =
            ViscosityData(ustrip.(T_data .|> K), [], ustrip.(ϱ), ustrip.(η), :unknown)

        @test all(data_eta_unit_pressure.T .== data_eta_ustrip_pressure.T)
        @test all(data_eta_unit_pressure.p .== data_eta_ustrip_pressure.p)
        @test all(data_eta_unit_pressure.Y .== data_eta_ustrip_pressure.Y)
        @test all(data_eta_unit_density.T .== data_eta_ustrip_density.T)
        @test all(data_eta_unit_density.Y .== data_eta_ustrip_density.Y)


        data_lamb_ustrip_pressure = ThermalConductivityData(
            ustrip.(T_data .|> K),
            ustrip.(P_data .|> Pa),
            [],
            ustrip.(λ),
            :unknown,
        )
        data_lamb_ustrip_density = ThermalConductivityData(
            ustrip.(T_data .|> K),
            [],
            ustrip.(ϱ),
            ustrip.(λ),
            :unknown,
        )

        @test all(data_lamb_unit_pressure.T .== data_lamb_ustrip_pressure.T)
        @test all(data_lamb_unit_pressure.p .== data_lamb_ustrip_pressure.p)
        @test all(data_lamb_unit_pressure.Y .== data_lamb_ustrip_pressure.Y)
        @test all(data_lamb_unit_density.T .== data_lamb_ustrip_density.T)
        @test all(data_lamb_unit_density.Y .== data_lamb_ustrip_density.Y)
    end

    ## Properties
    eos_model = PCSAFT(fluid)

    model_eta = FrameworkModel(eos_model, [data_eta_unit_density])
    model_lamb = FrameworkModel(eos_model, [data_lamb_unit_density])

    @testset "Viscosity" begin
        for (i, (p, t)) ∈ enumerate(PT)
            @test viscosity(model_eta, p, t) ≈ η[i] rtol = 1E-1
        end
    end

    @testset "Thermal Conductivity" begin
        for (i, (p, t)) ∈ enumerate(PT)
            @test thermal_conductivity(model_lamb, p, t) ≈ λ[i] rtol = 1E-1
        end
    end

    @testset "Unit Conversions" begin

        T = 25.0°C
        P = 1.0bar

        η_calc = viscosity(model_eta, P, T, output_unit = cP)
        λ_calc = thermal_conductivity(model_lamb, P, T, output_unit = btu / (hr * ft * Ra))

        @test η_calc ≈ 552.57μPa * s rtol = 1E-1
        @test λ_calc ≈ 130.35mW / m / K rtol = 1E-1

        P_psi = P |> psi
        T_Ra = T |> Ra

        η_calc2 = viscosity(model_eta, P_psi, T_Ra)
        λ_calc2 = thermal_conductivity(model_lamb, P_psi, T_Ra)

        @test η_calc2 ≈ η_calc rtol = 1E-13
        @test λ_calc2 ≈ λ_calc rtol = 1E-13
    end


    @testset "Transport Property Constructor" begin

        data_transp_p_eta = EntropyScaling.TransportPropertyData(T_data, P_data, η)
        data_transp_p_lamb = EntropyScaling.TransportPropertyData(T_data, P_data, λ)

        data_transp_ϱ_eta = EntropyScaling.TransportPropertyData(T_data, ϱ, η)
        data_transp_ϱ_lamb = EntropyScaling.TransportPropertyData(T_data, ϱ, λ)

        @test data_transp_p_eta.prop isa EntropyScaling.Viscosity
        @test data_transp_p_lamb.prop isa EntropyScaling.ThermalConductivity
        @test data_transp_ϱ_eta.prop isa EntropyScaling.Viscosity
        @test data_transp_ϱ_lamb.prop isa EntropyScaling.ThermalConductivity

        @test all(data_transp_p_eta.T .== data_eta_unit_pressure.T)
        @test all(data_transp_p_eta.p .== data_eta_unit_pressure.p)
        @test all(data_transp_p_eta.Y .== data_eta_unit_pressure.Y)
        @test all(data_transp_p_lamb.T .== data_lamb_unit_pressure.T)
        @test all(data_transp_p_lamb.p .== data_lamb_unit_pressure.p)
        @test all(data_transp_p_lamb.Y .== data_lamb_unit_pressure.Y)

        @test all(data_transp_ϱ_eta.T .== data_eta_unit_density.T)
        @test all(data_transp_ϱ_eta.Y .== data_eta_unit_density.Y)
        @test all(data_transp_ϱ_lamb.T .== data_lamb_unit_density.T)
        @test all(data_transp_ϱ_lamb.Y .== data_lamb_unit_density.Y)

    end

end
