
"""
    viscosity_CE(T, Mw, σ, ε)

Chapman-Enskog viscosity for the zero-density limit.
"""
function viscosity_CE(T, Mw, σ, ε)
    return 5/16 * √(Mw*kB*T/NA/π) / (σ^2*Ω_22(T*kB/ε))
end

"""
    viscosity_CE_plus(eos, T, σ, ε)

Scaled Chapman-Enskog viscosity for the zero-density limit.
"""
function viscosity_CE_plus(eos, T, σ, ε)
    dBdT = second_virial_coefficient_dT(eos,T)/NA
    B = second_virial_coefficient(eos,T)/NA
    return 5/16/√π / (σ^2*Ω_22(T*kB/ε)) * (T*dBdT+B)^(2/3)
end

"""
    thermal_conductivity_CE(T, Mw, σ, ε)

Chapman-Enskog thermal conductivity for the zero-density limit.
"""
function thermal_conductivity_CE(T, Mw, σ, ε)
    return 75/64 * kB * √(R*T/Mw/π) / (σ^2*Ω_22(T*kB/ε))
end

"""
    thermal_conductivity_CE_plus(eos, T, σ, ε)

Scaled Chapman-Enskog thermal conductivity for the zero-density limit.
"""
function thermal_conductivity_CE_plus(eos, T, σ, ε)
    dBdT = second_virial_coefficient_dT(eos,T)/NA
    B = second_virial_coefficient(eos,T)/NA
    return 75/64/√π / (σ^2*Ω_22(T*kB/ε)) * (T*dBdT+B)^(2/3)
end

"""
    self_diffusion_coefficient_CE(T, Mw, σ, ε)

Chapman-Enskog diffusion coefficient for the zero-density limit.
"""
function diffusion_coefficient_CE(T, Mw, σ, ε)
    return 3/8 * √(Mw*kB*T/NA/π) / (σ^2*Ω_11(T*kB/ε))
end

"""
    diffusion_coefficient_CE_plus(eos, T, σ, ε)

Scaled Chapman-Enskog diffusion coefficient for the zero-density limit.
"""
function diffusion_coefficient_CE_plus(eos, T, σ, ε)
    dBdT = second_virial_coefficient_dT(eos,T)/NA
    B = second_virial_coefficient(eos,T)/NA
    return 3/8/√π / (σ^2*Ω_11(T*kB/ε)) * (T*dBdT+B)^(2/3)
end

property_CE(prop::Viscosity, T, Mw, σ, ε) = viscosity_CE(T, Mw, σ, ε)
property_CE_plus(prop::Viscosity, eos, T, σ, ε) = viscosity_CE_plus(eos, T, σ, ε)
property_CE(prop::ThermalConductivity, T, Mw, σ, ε) = thermal_conductivity_CE(T, Mw, σ, ε)
property_CE_plus(prop::ThermalConductivity, eos, T, σ, ε) = thermal_conductivity_CE_plus(eos, T, σ, ε)
property_CE(prop::DiffusionCoefficient, T, Mw, σ, ε) = diffusion_coefficient_CE(T, Mw, σ, ε)
property_CE_plus(prop::DiffusionCoefficient, eos, T, σ, ε) = diffusion_coefficient_CE_plus(eos, T, σ, ε)

"""
    Ω_22(T_red)

Collision integrals from [Kim and Monroe (2014)](https://www.doi.org/10.1016/j.jcp.2014.05.018).
"""
function Ω_22(T_red)
    A = -0.92032979
    BCi = [ [   2.3508044,      1.6330213       ],
            [   0.50110649,     -6.9795156e-1   ],
            [   -4.7193769e-1,  1.6096572e-1,   ],
            [   1.5806367e-1,   -2.2109440e-2   ],
            [   -2.6367184e-2,  1.7031434e-3    ],
            [   1.8120118e-3,   -0.56699986e-4    ]]
    return A .+ sum([BCi[k][1] ./ T_red.^k .+ BCi[k][2] .* log.(T_red).^k for k in 1:6])
end

"""
    Ω_11(T_red)

Collision integrals from [Kim and Monroe (2014)](https://www.doi.org/10.1016/j.jcp.2014.05.018).
"""
function Ω_11(T_red)
    A = -1.1036729
    BCi = [ [   2.6431984,      1.6690746       ],
            [   0.0060432255,   -6.9145890e-1   ],
            [   -1.5158773e-1,  1.5502132e-1,   ],
            [   0.54237938e-1,  -2.0642189e-2   ],
            [   -0.90468682e-2, 1.5402077e-3    ],
            [   0.61742007e-3,  -0.49729535e-4  ]]
    return A .+ sum([BCi[k][1] ./ T_red.^k .+ BCi[k][2] .* log.(T_red).^k for k in 1:6])
end

"""
    correspondence_principle(Tc, pc)

Calculate the LJ parameters from the critical temperature and pressure.
"""
function correspondence_principle(Tc, pc)
    Tc_LJ = 1.321
    pc_LJ = 0.129

    ε = kB*Tc/Tc_LJ
    σ = (pc_LJ/pc*ε).^(1/3)
    return σ, ε
end