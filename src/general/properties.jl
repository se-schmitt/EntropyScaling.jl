export viscosity, thermal_conductivity, self_diffusion_coefficient
export MS_diffusion_coefficient, fick_diffusion_coefficient
export inf_diffusion_coefficient

const PROPERTY_FUNCTIONS = [
    :viscosity => :AbstractViscosity, 
    :thermal_conductivity => :AbstractThermalConductivity, 
    :self_diffusion_coefficient => :SelfDiffusionCoefficient,
    :MS_diffusion_coefficient => :MaxwellStefanDiffusionCoefficient, 
    :fick_diffusion_coefficient => :FickDiffusionCoefficient,
    :inf_diffusion_coefficient => :InfDiffusionCoefficient,
]

"""
    viscosity(model::EntropyScalingModel, p, T, z=[1.]; phase=:unknown)

Viscosity `η(p,T,x)` (`[η] = Pa s`).
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
    N = length(model)
    param = model[InfDiffusionCoefficient()]
    
    Ðᵢⱼ = zero(MSDiffusionMatrix, N)
    for i in 1:N, j in i+1:N
        #TODO extend to multicomponent mixtures
        s = entropy_conf(model.eos, ϱ, T, z)
        sˢ = scaling_variable(param, s, z)
        Dˢ = scaling_model(param, sˢ, z)
        Ðᵢⱼ[i,j] = scaling(param, model.eos, Dˢ, T, ϱ, s, z; inv=true)
    end
    return Ðᵢⱼ
end

"""
    fick_diffusion_coefficient(model::EntropyScalingModel, p, T, z; phase=:unknown)

Fickian diffusion coefficient `D(p,T,x)` (`[D] = m² s⁻¹`).
"""
fick_diffusion_coefficient

function fick_diffusion_coefficient(model::AbstractEntropyScalingModel, p, T, z; phase=:unknown)
    ϱ = molar_density(model.eos, p, T, z; phase=phase)
    return ϱT_fick_diffusion_coefficient(model, ϱ, T, z)
end

function ϱT_fick_diffusion_coefficient(model::AbstractEntropyScalingModel, ϱ, T, z)
    N = length(model)
    _rng = 1:N-1
    x = z ./ sum(z)
    Ð = ϱT_MS_diffusion_coefficient(model, ϱ, T, z)     #TODO multicomponent MS diffusion coefficient
    _Ð = inv.(Ð)
    setindex!.(Ref(_Ð), 0, 1:N, 1:N)
    Γ = thermodynamic_factor(model.eos, ϱ, T, z)
    B = [
        i == j ? 
        x[i]*_Ð[i,N] + sum(i == k ? 0 : x[k]*_Ð[i,k] for k in 1:N) :
        -x[i] * (_Ð[i,j] - _Ð[i,N])
    for i in _rng, j in _rng]
    D = inv(B) * Γ

    return D 
end

"""
    inf_diffusion_coefficient(model::EntropyScalingModel, p, T, z; phase=:unknown, solute=nothing, solvent=nothing)

Returns all diffusion coefficients at infinite dilution of the system (if parameters are available):

    [D₁  D₁₂ ⋯ D₁ₙ;
     D₂₁ D₂  ⋯ D₂ₙ;
     ⋮    ⋮   ⋱ ⋮  
     Dₙ₁ Dₙ₂ ⋯ Dₙ]

Dᵢⱼ is the diffusion coefficient of solute i at infinite dilution in solvent j.
If `solute` or `solvent` is specified, returns only the infinite diffusion coefficients in this component (one row or column of the matrix).
If both `solute` and `solvent` are specified, a scalar value is returned.
"""
inf_diffusion_coefficient

function inf_diffusion_coefficient(model::AbstractTransportPropertyModel, p, T; 
    phase=:unknown, solute=nothing, solvent=nothing
)
    
    TYPE = promote_type(typeof(p), typeof(T))
    N = length(model)
    idx_solute = isnothing(solute) ? (1:N) : match_comp(solute, model.components)
    idx_solvent = isnothing(solvent) ? (1:N) : match_comp(solvent, model.components)
   
    if all(length.([idx_solute,idx_solvent]) .== 1)
        idx_i, idx_j = only(idx_solute), only(idx_solvent)
        Dij = _inf_diffusion_coefficient(model, p, T, (idx_i, idx_j); phase)
    else
        Dij = zeros(TYPE, length(idx_solute), length(idx_solvent))
        for (i,idx_i) in enumerate(idx_solute), (j,idx_j) in enumerate(idx_solvent)
            if i != j
                Dij[i,j] = _inf_diffusion_coefficient(model, p, T, (idx_i, idx_j); phase)
            end
        end
    end
    return Dij
end

match_comp(comp::AbstractString, components) = findall(comp .== components)
match_comp(comp::Int, components) = [comp]

function _inf_diffusion_coefficient(model, p, T, (idx_i, idx_j); phase)
end