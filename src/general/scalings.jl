
"""
    rosenfeld_scaling(param::BaseParam{Viscosity},η,T,ϱN;inv=false)
    rosenfeld_scaling(param::BaseParam{ThermalConductivity},λ,T,ϱN;inv=false)
    rosenfeld_scaling(param::BaseParam{DiffusionCoefficient},D,T,ϱN;inv=false)

Rosenfeld scaling for transport properties.

## References

"""
rosenfeld_scaling

function rosenfeld_scaling(param::BaseParam{<:AbstractViscosity}, η, T, ϱ, z = Z1; inv=false)
    k = !inv ? 1 : -1
    Mw = _dot(param.Mw,z)
    return η / ((ϱ*NA)^(2/3) * sqrt(Mw / NA * kB * T))^k
end

function rosenfeld_scaling(param::BaseParam{<:AbstractThermalConductivity}, λ, T, ϱ, z = Z1; inv=false)
    k = !inv ? 1 : -1
    Mw = _dot(param.Mw,z)
    return λ / ((ϱ*NA)^(2/3) * kB)^k * sqrt(Mw / (T * kB * NA))^k
end

function rosenfeld_scaling(param::BaseParam{<:DiffusionCoefficient}, D, T, ϱ, z = Z1; inv=false)
    k = !inv ? 1 : -1
    Mw = _dot(param.Mw,z)
    return D * sqrt(Mw / (NA * kB * T))^k * (ϱ*NA)^(k/3)
end

"""
    plus_scaling(param::BaseParam,Y,T,ϱN,s;inv=false)

Plus scaling for transport properties.

## References   
"""
plus_scaling

function plus_scaling(param::BaseParam, Y, T, ϱ, s, z = Z1; inv=false)
    k = !inv ? 1 : -1
    return rosenfeld_scaling(param, Y, T, ϱ, z; inv=inv) * (-s/R)^(k*2/3)
end