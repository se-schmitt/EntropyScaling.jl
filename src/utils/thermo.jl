second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> CL.second_virial_coefficient(eos, xT, z), T)

get_m(eos::Any) = fun_na_error("get_m",typeof(eos))
get_components(eos::Any) = fun_na_error("get_components",typeof(eos))

# Bulk properties
pressure(eos::CL.EoSModel, ϱ, T, z=Z1) = CL.pressure(eos, 1.0/ϱ, T, z)

molar_density(eos::CL.EoSModel, p, T, z=Z1; phase=:unknown, ϱ0=nothing) =
    CL.molar_density(eos, p, T, z; phase, vol0=isnothing(ϱ0) ? nothing : 1/ϱ0)

molar_density(eos::CL.EoSVectorParam, p, T, z=Z1; phase=:unknown, ϱ0=nothing) =
    molar_density(eos.model, p, T, z; phase, ϱ0)

# note: each EoSModel in Clapeyron can set its own gas constant; divide by that and multiply
# by EntropyScaling's R so calculations are consistent between models
entropy_conf(eos::CL.EoSModel, ϱ, T, z=Z1) =
    CL.VT_entropy_res(eos, 1.0/ϱ, T, z) * R / CL.Rgas(eos)
entropy_conf(eos::CL.EoSVectorParam, ϱ, T, z=Z1) = entropy_conf(eos.model, ϱ, T, z)

second_virial_coefficient(eos::CL.EoSModel, T, z=Z1) = CL.second_virial_coefficient(eos, T, z)
second_virial_coefficient(eos::CL.EoSVectorParam, T, z=Z1) = second_virial_coefficient(eos.model, T, z)

isobaric_heat_capacity(eos::CL.EoSModel, ϱ, T, z=Z1) = CL.VT_isobaric_heat_capacity(eos, 1.0/ϱ, T, z)
isobaric_heat_capacity(eos::CL.EoSVectorParam, ϱ, T, z=Z1) = isobaric_heat_capacity(eos.model, ϱ, T, z)

isochoric_heat_capacity(eos::CL.EoSModel, ϱ, T, z=Z1) = CL.VT_isochoric_heat_capacity(eos, 1.0/ϱ, T, z)
isochoric_heat_capacity(eos::CL.EoSVectorParam, ϱ, T, z=Z1) = isochoric_heat_capacity(eos.model, ϱ, T, z)

thermodynamic_factor(eos::CL.EoSModel, ϱ, T, z=Z1) = CL.VT_thermodynamic_factor(eos, 1.0/ϱ, T, z)
thermodynamic_factor(eos::CL.EoSVectorParam, ϱ, T, z=Z1) = CL.VT_thermodynamic_factor(eos.model, 1.0/ϱ, T, z)

# Critical properties
function crit_pure(eos::CL.EoSModel)
    Tc, pc, vc = CL.crit_pure(eos)
    return Tc, pc, 1/vc
end
function crit_mix(eos::CL.EoSModel, z)
    Tc, pc, vc = CL.crit_mix(eos, z)
    return Tc, pc, 1/vc
end

get_m(eos::CL.EoSModel) = :segment in fieldnames(typeof(eos.params)) ? eos.params.segment.values : ones(length(eos.name))
get_m(eos::CL.EoSVectorParam) = get_m(eos.model)
get_m(eos::CL.SingleFluid) = SA1
get_m(eos::CL.SAFTgammaMieModel) = get_m(eos.vr_model)

get_components(eos::CL.EoSModel) = eos.components

# CL.init_model extensions (used by JutulDarcy integration)
function CL.init_model(model::AESM, components, userlocations=String[], verbose=false, reference_state=nothing)
    return model
end
function CL.init_model(::Type{𝕄}, components, userlocations=String[], verbose=false, reference_state=nothing) where {𝕄 <: AESM}
    return 𝕄(components)
end
