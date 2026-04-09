const SA1 = CL.SA[1.0] #TODO remove

# Function to throw error
fun_na_error(fname,eostype) = error("Function `$fname` needs to be defined for eos type `$eostype`.")

# Bulk properties
"""
    pressure(eos, ϱ, T, z=[1.])

Pressure `p` (`[p] = Pa`).
"""
pressure(eos::Any, ϱ, T, z=Z1) = fun_na_error("pressure",typeof(eos))

"""
    molar_density(eos, p, T, z=[1.]; phase=:unknown)

Molar density `ϱ` (`[ϱ] = mol/m³`).
"""
molar_density(eos::Any, p, T, z=Z1; phase=:unknown, ϱ0=nothing) = fun_na_error("molar_density",typeof(eos))

"""
    entropy_conf(eos, ϱ, T, z=[1.])

Configurational (or residual) entropy `s_conf` (`[s_conf] = J/(mol K)`).
"""
entropy_conf(eos::Any, ϱ, T, z=Z1) = fun_na_error("entropy_conf",typeof(eos))

"""
    second_virial_coefficient(eos, T, z=[1.])

Second virial coefficient `B` (`[B] = m³/mol`).
"""
second_virial_coefficient(eos::Any, T, z=Z1) = fun_na_error("second_virial_coefficient",typeof(eos))

"""
    second_virial_coefficient_dT(eos, T, z=[1.])

Temperature derivative of the second virial coefficient `dBdT` (`[dBdT] = m³/(mol K)`).
"""
second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> second_virial_coefficient(eos, xT, z), T)

"""
    isobaric_heat_capacity(eos, ϱ, T, z=[1.])

Isobaric heat capacity `cₚ` (`[cₚ] = J/(mol K)`).
"""
isobaric_heat_capacity(eos::Any, ϱ, T, z=Z1) = fun_na_error("isobaric_heat_capacity",typeof(eos))

"""
    isochoric_heat_capacity(eos, ϱ, T, z=[1.])

Isochoric heat capacity `cᵥ` (`[cᵥ] = J/(mol K)`).
"""
isochoric_heat_capacity(eos::Any, ϱ, T, z=Z1) = fun_na_error("isochoric_heat_capacity",typeof(eos))

"""
    thermodynamic_factor(eos, ϱ, T, z)

Thermodynamic factor `Γᵢⱼ`.
"""
thermodynamic_factor(eos::Any, ϱ, T, z) = fun_na_error("thermdynamic_factor",typeof(eos))

# Critical properties
"""
    crit_pure(eos)

Critical point of pure component (`Tc`, `pc`, `ϱc`) (`[Tc] = K`, `[pc] = Pa`, `[ϱc] = mol/m³`).
"""
crit_pure(eos::Any) = fun_na_error("crit_pure",typeof(eos))

"""
    crit_mix(eos, x)

Critical point of a mixture at given composition (`Tc`, `pc`, `ϱc`) (`[Tc] = K`, `[pc] = Pa`, `[ϱc] = mol/m³`).
"""
crit_mix(eos::Any, z) = fun_na_error("crit_mix",typeof(eos))

# Utility functions
"""
    split_model(eos)

Split the model in pure component models.
"""
split_model(eos::Any) = fun_na_error("split_model",typeof(eos))


"""
    get_m(eos)

Segment parameter of SAFT-type EOS.
"""
get_m(eos::Any) = fun_na_error("get_m",typeof(eos))

"""
    get_components(eos)

Component names of the system.
"""
get_components(eos::Any) = fun_na_error("get_components",typeof(eos))

# Bulk properties
pressure(eos::CL.EoSModel, ϱ, T, z=SA1) = CL.pressure(eos, 1.0/ϱ, T, z)

molar_density(eos::CL.EoSModel, p, T, z=SA1; phase=:unknown, ϱ0=nothing) =
    CL.molar_density(eos, p, T, z; phase, vol0=isnothing(ϱ0) ? nothing : 1/ϱ0)

molar_density(eos::CL.EoSVectorParam, p, T, z=SA1; phase=:unknown, ϱ0=nothing) =
    molar_density(eos.model, p, T, z; phase, ϱ0)

# note: each EoSModel in Clapeyron can set its own gas constant; divide by that and multiply
# by EntropyScaling's R so calculations are consistent between models
entropy_conf(eos::CL.EoSModel, ϱ, T, z=SA1) =
    CL.VT_entropy_res(eos, 1.0/ϱ, T, z) * R / CL.Rgas(eos)
entropy_conf(eos::CL.EoSVectorParam, ϱ, T, z=SA1) = entropy_conf(eos.model, ϱ, T, z)

second_virial_coefficient(eos::CL.EoSModel, T, z=SA1) = CL.second_virial_coefficient(eos, T, z)
second_virial_coefficient(eos::CL.EoSVectorParam, T, z=SA1) = second_virial_coefficient(eos.model, T, z)

isobaric_heat_capacity(eos::CL.EoSModel, ϱ, T, z=SA1) = CL.VT_isobaric_heat_capacity(eos, 1.0/ϱ, T, z)
isobaric_heat_capacity(eos::CL.EoSVectorParam, ϱ, T, z=SA1) = isobaric_heat_capacity(eos.model, ϱ, T, z)

isochoric_heat_capacity(eos::CL.EoSModel, ϱ, T, z=SA1) = CL.VT_isochoric_heat_capacity(eos, 1.0/ϱ, T, z)
isochoric_heat_capacity(eos::CL.EoSVectorParam, ϱ, T, z=SA1) = isochoric_heat_capacity(eos.model, ϱ, T, z)

thermodynamic_factor(eos::CL.EoSModel, ϱ, T, z=SA1) = CL.VT_thermodynamic_factor(eos, 1.0/ϱ, T, z)
thermodynamic_factor(eos::CL.EoSVectorParam, ϱ, T, z=SA1) = CL.VT_thermodynamic_factor(eos.model, 1.0/ϱ, T, z)

# Critical properties
function crit_pure(eos::CL.EoSModel)
    Tc, pc, vc = CL.crit_pure(eos)
    return Tc, pc, 1/vc
end
function crit_mix(eos::CL.EoSModel, z)
    Tc, pc, vc = CL.crit_mix(eos, z)
    return Tc, pc, 1/vc
end

# Utility functions
split_model(eos::CL.EoSModel) = CL.split_model(eos)
split_model(eos::CL.SingleFluid) = [eos]
split_model(eos::CL.MultiFluid) = eos.pures
split_model(eos::CL.EoSVectorParam) = eos.pure

get_Mw(eos::CL.EoSModel) = CL.mw(eos)
get_Mw(eos::CL.EoSVectorParam) = get_Mw(eos.model)

get_m(eos::CL.EoSModel) = :segment in fieldnames(typeof(eos.params)) ? eos.params.segment.values : ones(length(eos.name))
get_m(eos::CL.EoSVectorParam) = get_m(eos.model)
get_m(eos::CL.SingleFluid) = SA1
get_m(eos::CL.SAFTgammaMieModel) = get_m(eos.vr_model)

get_components(eos::CL.EoSModel) = eos.components

_eos_cache(eos) = eos
_eos_cache(eos::CL.EoSModel) = CL.EoSVectorParam(eos)
_eos_cache(eos::CL.EoSVectorParam) = eos
_eos_cache(eos::CL.MultiFluid) = eos

# CL.init_model extensions (used by JutulDarcy integration)
function CL.init_model(model::AESM, components, userlocations=String[], verbose=false, reference_state=nothing)
    return model
end
function CL.init_model(::Type{𝕄}, components, userlocations=String[], verbose=false, reference_state=nothing) where {𝕄 <: AESM}
    return 𝕄(components)
end
