using Test
using EntropyScaling
using Clapeyron


model_but = PCSAFT("n-butane")

@testset "fit" begin
    # Viscosity
    exp_η = [
        298.15000000 1.36392960 0.00000738;
        323.14000000 1.27393424 0.00000799;
        344.26000000 2.74292056 0.00000850;
        360.92000000 2.72966139 0.00000914;
        410.92000000 1.97974404 0.00001018;
        494.26000000 234.50438089 0.00002400;
        448.16000000 321.10342549 0.00003694;
        410.93000000 426.61955224 0.00006300;
        360.93000000 509.95555500 0.00010050;
        327.59000000 543.08771615 0.00012650;
        413.00000000 560.34725853 0.00013460;
        410.93000000 570.85970426 0.00014550;
        250.00000000 651.91452465 0.00033730;
        140.00000000 734.65968124 0.00201250;
        140.00000000 744.80258283 0.00238460;
    ]
    (α_η, ηˢ, s_η) = fit_entropy_scaling(model_but, exp_η[:,1], exp_η[:,2], exp_η[:,3], "vis")
    @test all(isapprox(α_η, [0.0,-9.490597,9.747817,-1.279090,0.3666153]; atol=1e-5))

    # Thermal conductivity
    exp_λ = [
        264.650000 0.026422 0.012470;
        260.010000 0.815645 0.013130;
        284.240000 0.769328 0.015360;
        291.420000 2.543079 0.015360;
        291.420000 2.593972 0.015360;
        291.920000 2.614573 0.015360;
        420.900000 366.248883 0.067940;
        333.570000 607.841347 0.122700;
        281.810000 631.047997 0.131000;
        230.950000 651.424554 0.140900;
        261.900000 671.272085 0.152600;
        262.550000 670.756747 0.152600;
        190.350000 686.052327 0.158100;
        211.160000 683.615938 0.158100;
        211.670000 683.182668 0.158100;
        170.070000 707.993986 0.167700;
        141.540000 733.494966 0.175300;
        170.120000 720.567767 0.175300;
        150.720000 756.096879 0.188900;
    ]
    (α_λ, λˢ, s_λ) = fit_entropy_scaling(model_but, exp_λ[:,1], exp_λ[:,2], exp_λ[:,3], "tcn"; i_fit=[1,1,1,1,1])
    @test all(isapprox(α_λ, [3.532197, 116.820497, -99.444059, 24.246002, 0.58713]; atol=1e-5))

    # Self-diffusion coefficient
    exp_D = [
        451.00 370.72 0.00000003;
        361.00 503.20 0.0000000134;
        451.00 628.38 0.00000000901;
        398.00 616.33 0.00000000811;
        361.00 635.74 0.00000000655;
        293.00 584.95 0.00000000623;
        361.00 667.37 0.00000000486;
        293.00 673.01 0.00000000345;
        303.00 716.76 0.00000000265;
        293.00 721.37 0.00000000226;
        235.00 731.63 0.00000000168;
        177.00 696.62 0.00000000125;
        177.00 725.40 9.51E-10;
        203.00 771.06 6.73E-10;
        150.00 745.15 4.68E-10;
        177.00 789.99 3.13E-10;
        150.00 765.47 2.95E-10;
    ]
    (α_D, Dˢ, s_D) = fit_entropy_scaling(model_but, exp_D[:,1], exp_D[:,2], exp_D[:,3], "dif"; i_fit=[0,0,0,1,1])
    @test all(isapprox(α_D, [0.0, 0.0, 0.0, -3.573438, -0.99224]; atol=1e-5))
end

model_tol = PCSAFT("toluene")
@testset "fit inf dil" begin
    # Infinite dilution diffusion coefficient
    exp_D_inf = [
        299.15 851.379445 0.00000000252;
        299.15 858.494238 0.00000000234;
        299.15 866.073199 0.00000000223;
        299.15 873.609644 0.00000000198;
        323.15 829.423562 0.00000000321;
        323.15 837.290948 0.00000000314;
        323.15 845.315219 0.00000000299;
        323.15 853.174715 0.0000000028;
        348.15 806.330239 0.00000000417;
        348.15 815.643854 0.00000000401;
        348.15 823.880471 0.0000000038;
        348.15 833.2534 0.00000000359;
    ]
    solute = Dict(:M=>0.08618, :Tc=>519.33427, :pc=>3.542718e6)
    (α_D_inf, D_infˢ, s_D_inf) = fit_entropy_scaling(model_but, exp_D_inf[:,1], exp_D_inf[:,2], exp_D_inf[:,3], "dif"; solute=solute, i_fit=[0,0,0,1,1])   
    @test all(isapprox(α_D_inf, [0.0, 0.0, 0.0, -1.153974, -0.5866215]; atol=1e-5))
end

model_mix1 = PCSAFT(["benzene","hexane"])
model_mix2 = PCPSAFT(["acetone", "chloroform"]; userlocations=(;dipole=[2.88,0.0],epsilon=[232.99,271.63],sigma=[3.2742,3.4709],Mw=[58.08,119.38],segment=[2.7447,2.5038]))
model_mix3 = PCSAFT(["toluene","heptane"])

α_η_but = [[0.,-14.165,13.97,-2.382,0.501]]
α_λ_but = [[3.962,98.222,-82.974,20.079,1.073]]
α_D_but = [[0.,0.,0.,-3.507,-0.997]]
α_λ_mix1 = [
    [6.492054425320112,0.0,0.0,1.9855528414772259,3.12643833453098],
    [9.75696539716926,0.0,0.0,1.6572105176108498,5.058427863917941]
]
α_D_mix2 = [
    [0.,0.,0.,1.1301081958e+01,-7.7638664176e+00],
    [0.,0.,0.,-3.0212807487e+00,-1.8500112748e+00]
]
α_Ð_mix3 = [
    [0.,0.,0.,1.1606624185e+01,-7.4861787912e+00],
    [0.,0.,0.,1.8829820233e+01,-1.3186835311e+01],
]
@testset "call" begin
    @test call_entropy_scaling(model_but,  [323.],   [ 602.4253281], "vis"; α=α_η_but)[1] ≈ 1.9203575445e-4 rtol=1e-3
    @test call_entropy_scaling(model_but,  [323.],   [ 602.4253281], "tcn"; α=α_λ_but)[1] ≈ 0.11980321582 rtol=1e-3
    @test call_entropy_scaling(model_but,  [323.],   [ 602.4253281], "dif"; α=α_D_but)[1] ≈ 6.6470158395e-9 rtol=1e-3
    @test call_entropy_scaling(model_mix1, [294.7],  [ 740.7881963], "tcn"; x=[.5 .5;], α=α_λ_mix1)[1] ≈ 0.1338505 rtol=1e-3
    @test call_entropy_scaling(model_mix2, [298.15], [1151.7720905], "selfdif"; x=[.5 .5;], α=α_D_mix2, difcomp=1)[1] ≈ 2.872998e-9 rtol=1e-3
    @test call_entropy_scaling(model_mix3, [308.15], [ 741.0969526], "mutdif"; x=[.5 .5;], α=α_Ð_mix3)[1] ≈ 3.3341968265e-9 rtol=1e-3
end