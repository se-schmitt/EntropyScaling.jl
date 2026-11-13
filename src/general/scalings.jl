
"""
    rosenfeld_scaling(param::BaseParam{Viscosity}, η, T, ϱ; inverse=false)
    rosenfeld_scaling(param::BaseParam{ThermalConductivity}, λ, T, ϱ; inverse=false)
    rosenfeld_scaling(param::BaseParam{DiffusionCoefficient}, D, T, ϱ; inverse=false)

Rosenfeld scaling for transport properties.

## References

"""
rosenfeld_scaling

function rosenfeld_scaling(param::BaseParam{<:AbstractViscosity}, η, T, ϱ, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    (; kB, NA) = param.units
    Mw = _dot(param.Mw,z)*param.units.Mw
    return η / ((ϱ*NA)^(2/3) * sqrt(Mw / NA * kB * T))^k
end

function rosenfeld_scaling(param::BaseParam{<:AbstractThermalConductivity}, λ, T, ϱ, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    (; kB, NA) = param.units
    Mw = _dot(param.Mw,z)*param.units.Mw
    return λ / ((ϱ*NA)^(2/3) * kB)^k * sqrt(Mw / (T * kB * NA))^k
end

function rosenfeld_scaling_mass_frac_mix(param::BaseParam{<:AbstractThermalConductivity}, λ, T, ϱ, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    
    Mw = collect(param.Mw.values)*param.units.Mw
    m = Mw / NA
    y = (z .* Mw) / sum(z .* Mw)
    mₘᵢₓ = sum(y .* sqrt.(m))^2

    return λ / ((ϱ*NA)^(2/3) * kB)^k * sqrt(mₘᵢₓ/(kB*T))^k
end

function rosenfeld_scaling(param::BaseParam{<:AbstractDiffusionCoefficient}, D, T, ϱ, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    (; kB, NA) = param.units
    Mw = _dot(param.Mw,z)*param.units.Mw
    return D * sqrt(Mw / (NA * kB * T))^k * (ϱ*NA)^(k/3)
end

"""
    plus_scaling(param::BaseParam, Y, T, ϱ, s; inverse=false)

Plus scaling for transport properties.

## References   
"""
plus_scaling

function plus_scaling(param::BaseParam, Y, T, ϱ, s, z = Z1; inverse=false, mass_frac_weighted_mixing=false)
    k = !inverse ? 1 : -1
    if !mass_frac_weighted_mixing
        return rosenfeld_scaling(param, Y, T, ϱ, z; inverse) * (-s/(param.units.R))^(k*2/3)
    else
        return rosenfeld_scaling_mass_frac_mix(param, Y, T, ϱ, z; inverse) * (-s/(param.units.R))^(k*2/3)
    end
end