
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

property_CE(prop::AbstractViscosity, T, Mw, σ, ε) = viscosity_CE(T, Mw, σ, ε)
property_CE_plus(prop::AbstractViscosity, eos, T, σ, ε) = viscosity_CE_plus(eos, T, σ, ε)
property_CE(prop::AbstractThermalConductivity, T, Mw, σ, ε) = thermal_conductivity_CE(T, Mw, σ, ε)
property_CE_plus(prop::AbstractThermalConductivity, eos, T, σ, ε) = thermal_conductivity_CE_plus(eos, T, σ, ε)
property_CE(prop::DiffusionCoefficient, T, Mw, σ, ε) = diffusion_coefficient_CE(T, Mw, σ, ε)
property_CE_plus(prop::DiffusionCoefficient, eos, T, σ, ε) = diffusion_coefficient_CE_plus(eos, T, σ, ε)

"""
    Ω_22(T_red)

Collision integrals from [Kim and Monroe (2014)](https://www.doi.org/10.1016/j.jcp.2014.05.018).
"""
function Ω_22(T_red)
    A = -0.92032979
    BCi = ( (   2.3508044,      1.6330213      ),
            (   0.50110649,     -6.9795156e-1  ),
            (   -4.7193769e-1,  1.6096572e-1,  ),
            (   1.5806367e-1,   -2.2109440e-2  ),
            (   -2.6367184e-2,  1.7031434e-3   ),
            (   1.8120118e-3,   -0.56699986e-4 ))
    Ω22 = zero(1.0*T_red) + A
    T_red_i = T_red
    lnT_red = log(T_red)
    lnT_red_i = lnT_red
    for k in 1:6
        BCi1,BCi2 = BCi[k]
        Ω22 += BCi1/T_red_i + BCi2*lnT_red_i
        T_red_i *= T_red
        lnT_red_i *= lnT_red
    end
    return Ω22
    #return A .+ sum([BCi[k][1] ./ T_red.^k .+ BCi[k][2] .* log.(T_red).^k for k in 1:6])
end

"""
    Ω_11(T_red)

Collision integrals from [Kim and Monroe (2014)](https://www.doi.org/10.1016/j.jcp.2014.05.018).
"""
function Ω_11(T_red)
    A = -1.1036729
    BCi = ( (   2.6431984,      1.6690746      ),
            (   0.0060432255,   -6.9145890e-1  ),
            (   -1.5158773e-1,  1.5502132e-1,  ),
            (   0.54237938e-1,  -2.0642189e-2  ),
            (   -0.90468682e-2, 1.5402077e-3   ),
            (   0.61742007e-3,  -0.49729535e-4 ))
    Ω11 = zero(1.0*T_red) + A
    T_red_i = T_red
    lnT_red = log(T_red)
    lnT_red_i = lnT_red
    for k in 1:6
        BCi1,BCi2 = BCi[k]
        Ω11 += BCi1/T_red_i + BCi2*lnT_red_i
        T_red_i *= T_red
        lnT_red_i *= lnT_red
    end
    return Ω11
    #return A .+ sum([BCi[k][1] ./ T_red.^k .+ BCi[k][2] .* log.(T_red).^k for k in 1:6])
end

"""
    correspondence_principle(Tc, pc)

Calculate the LJ parameters from the critical temperature and pressure.
"""
function correspondence_principle(Tc, pc)
    Tc_LJ = 1.321
    pc_LJ = 0.129

    ε = kB*Tc/Tc_LJ
    σ = cbrt(pc_LJ/pc*ε)
    return σ, ε
end

# Mixing rule from Wilke (1950) [DOI: 10.1063/1.1747673] and Mason and Saxena (1958) [DOI: 10.1063/1.1724352]
"""
    mix_CE(param::BaseParam{Viscosity}, x)
    mix_CE(param::BaseParam{ThermalConductivity}, x)

Mixing rule for Chapman-Enskog transport properties by Wilke (1950) for Viscosity and by 
Mason and Saxena (1958) for ThermalConductivity.

## References
(1) Wilke, C. R. A Viscosity Equation for Gas Mixtures. The Journal of Chemical Physics 
1950, 18 (4), 517–519. https://doi.org/10.1063/1.1747673.
(2) Mason, E. A.; Saxena, S. C. Approximate Formula for the Thermal Conductivity of Gas 
Mixtures. The Physics of Fluids 1958, 1 (5), 361–369. https://doi.org/10.1063/1.1724352.
"""
function mix_CE(param::BaseParam{P}, Y, x) where {P <: Union{AbstractViscosity, AbstractThermalConductivity}}
    Y₀_mix = zero(Base.promote_eltype(param.T_range,Y,x))
    zero(Base.promote_eltype(param.T_range,Y,x))
    enum_M = enumerate(param.Mw)
    for (i,Mi) in enum_M
        xΦ = zero(Y₀_mix)
        for (j,Mj) in enum_M
            xΦ += x[j] * (1+√(Y[i]/Y[j])*√√(Mj/Mi))^2 / √(8*(1+Mi/Mj))
        end
        #xΦ = sum([x[j] * (1+√(Y[i]/Y[j])*√√(Mj/Mi))^2 / √(8*(1+Mi/Mj)) for (j,Mj) in enum_M])
        Y₀_mix += x[i]*Y[i]/xΦ
    end
    return Y₀_mix
end

"""
    mix_CE(param::BaseParam{DiffusionCoefficient}, x)

Mixing rule for Chapman-Enskog diffusion coefficient by Miller and Carman (1961).

## References
(1) Miller, L.; Carman, P. C. Self-Diffusion in Mixtures. Part 4. -- Comparison of Theory 
and Experiment for Certain Gas Mixtures. Trans. Faraday Soc. 1961, 57 (0), 2143–2150. 
    https://doi.org/10.1039/TF9615702143.
"""
function mix_CE(param::BaseParam{P}, Y, x) where {P <: DiffusionCoefficient}
    return 1.0 / sum(x[i] / Y[i] for i in eachindex(Y))
end

calc_M_CE(Mw) = 2.0/sum(inv,Mw)