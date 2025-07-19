"""
    viscosity(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Dynamic viscosity `η(p,T,x)` (`[η] = Pa s`).
"""
viscosity

function viscosity(model::AbstractEntropyScalingModel, p, T, z=Z1; phase=:unknown)
    ϱ = molar_density(model.eos, p, T, z; phase=phase)
    return ϱT_viscosity(model, ϱ, T, z)
end

function ϱT_viscosity(model::AbstractEntropyScalingModel, ϱ, T, z::AbstractVector=Z1)
    param = model[Viscosity()]
    s = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)
    ηˢ = scaling_model(param, sˢ, z)
    return scaling(param, model.eos, ηˢ, T, ϱ, s, z; inv=true)
end


"""
    thermal_conductivity(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Thermal conductivity `λ(p,T,x)` (`[λ] = W m⁻¹ K⁻¹`).
"""
thermal_conductivity

function thermal_conductivity(model::AbstractEntropyScalingModel, p, T, z=Z1; phase=:unknown)
    ϱ = molar_density(model.eos, p, T, z; phase=phase)
    return ϱT_thermal_conductivity(model, ϱ, T, z)
end

function ϱT_thermal_conductivity(model::AbstractEntropyScalingModel, ϱ, T, z::AbstractVector=Z1)
    param = model[ThermalConductivity()]
    s = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)
    λˢ = scaling_model(param, sˢ, z)
    return scaling(param, model.eos, λˢ, T, ϱ, s, z; inv=true)
end

"""
    self_diffusion_coefficient(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Self-diffusion coefficient `D(p,T,x)` (`[D] = m² s⁻¹`).
"""
self_diffusion_coefficient

function self_diffusion_coefficient(model::AbstractEntropyScalingModel, p, T, z=Z1; phase=:unknown)
    ϱ = molar_density(model.eos, p, T, z; phase=phase)
    if length(model) == 1
        return ϱT_self_diffusion_coefficient(model, ϱ, T)
    else
        return ϱT_self_diffusion_coefficient(model, ϱ, T, z)
    end
end

function ϱT_self_diffusion_coefficient(model::AbstractEntropyScalingModel, ϱ, T)
    param = model[SelfDiffusionCoefficient()]
    s = entropy_conf(model.eos, ϱ, T)
    sˢ = scaling_variable(param, s)
    Dˢ = scaling_model(param, sˢ)
    return scaling(param, model.eos, Dˢ, T, ϱ, s; inv=true)
end

"""
    MS_diffusion_coefficient(model::EntropyScalingModel, p, T, z; phase=:unknown)

Maxwell-Stefan diffusion coefficient `Ð(p,T,x)` (`[Ð] = m² s⁻¹`).
"""
MS_diffusion_coefficient

function MS_diffusion_coefficient(model::AbstractEntropyScalingModel, p, T, z; phase=:unknown)
    ϱ = molar_density(model.eos, p, T, z; phase=phase)
    return ϱT_MS_diffusion_coefficient(model, ϱ, T, z)
end

function ϱT_MS_diffusion_coefficient(model::AbstractEntropyScalingModel, ϱ, T, z)
    param = model[InfDiffusionCoefficient()]
    s = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)
    Dˢ = scaling_model(param, sˢ, z)
    return scaling(param, model.eos, Dˢ, T, ϱ, s, z; inv=true)
end