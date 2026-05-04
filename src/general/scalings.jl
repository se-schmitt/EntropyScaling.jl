
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
    Mw = _dot(param.Mw,z)*1e-3
    return η / ((ϱ*NA)^(2/3) * sqrt(Mw / NA * kB * T))^k
end

function rosenfeld_scaling(param::BaseParam{<:AbstractThermalConductivity}, λ, T, ϱ, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    Mw = _dot(param.Mw,z)*1e-3
    return λ / ((ϱ*NA)^(2/3) * kB)^k * sqrt(Mw / (T * kB * NA))^k
end

function rosenfeld_scaling(param::BaseParam{<:AbstractDiffusionCoefficient}, D, T, ϱ, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    Mw = _dot(param.Mw,z)*1e-3
    return D * sqrt(Mw / (NA * kB * T))^k * (ϱ*NA)^(k/3)
end

"""
    plus_scaling(param::BaseParam, Y, T, ϱ, s; inverse=false)

Plus scaling for transport properties.

## References   
"""
plus_scaling

function plus_scaling(param::BaseParam, Y, T, ϱ, s, z = Z1; inverse=false)
    k = !inverse ? 1 : -1
    return rosenfeld_scaling(param, Y, T, ϱ, z; inverse) * (-s/R)^(k*2/3)
end