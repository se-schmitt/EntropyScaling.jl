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

Molar density `ϱ` (`[ϱ] = mol m⁻³`).
"""
molar_density(eos::Any, p, T, z=Z1; phase=:unknown, ϱ0=nothing) = fun_na_error("molar_density",typeof(eos))

"""
    entropy_conf(eos, ϱ, T, z=[1.])

Configurational (or residual) entropy `s_conf` (`[s_conf] = J (mol K)⁻¹`).
"""
entropy_conf(eos::Any, ϱ, T, z=Z1) = fun_na_error("entropy_conf",typeof(eos))

"""
    second_virial_coefficient(eos, T, z=[1.])

Second virial coefficient `B` (`[B] = m³ mol⁻¹`).
"""
second_virial_coefficient(eos::Any, T, z=Z1) = fun_na_error("second_virial_coefficient",typeof(eos))

"""
    second_virial_coefficient_dT(eos, T, z=[1.])

Temperature derivative of the second virial coefficient `dBdT` (`[dBdT] = m³ (mol K)⁻¹`).
"""
second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> second_virial_coefficient(eos, xT, z), T)

"""
    isobaric_heat_capacity(eos, ϱ, T, z=[1.])

Isobaric heat capacity `cₚ` (`[cₚ] = J (mol K)⁻¹`).
"""
isobaric_heat_capacity(eos::Any, ϱ, T, z=Z1) = fun_na_error("isobaric_heat_capacity",typeof(eos))

"""
    isochoric_heat_capacity(eos, ϱ, T, z=[1.])

Isochoric heat capacity `cₚ` (`[cₚ] = J (mol K)⁻¹`).
"""
isochoric_heat_capacity(eos::Any, ϱ, T, z=Z1) = fun_na_error("isochoric_heat_capacity",typeof(eos))

# Critical properties
"""
    crit_pure(eos)

Critical point of pure component (`Tc`, `pc`, `ϱc`) (`[Tc] = K`, `[pc] = Pa`, `[ϱc] = mol m⁻³`).
"""
crit_pure(eos::Any) = fun_na_error("crit_pure",typeof(eos))

"""
    crit_mix(eos, x)

Critical point of a mixture at given composition (`Tc`, `pc`, `ϱc`) (`[Tc] = K`, `[pc] = Pa`, `[ϱc] = mol m⁻³`).
"""
crit_mix(eos::Any, z) = fun_na_error("crit_mix",typeof(eos))

# Utility functions
"""
    split_model(eos)

Split the model in pure component models.
"""
split_model(eos::Any) = fun_na_error("split_model",typeof(eos))

"""
    get_Mw(eos)

Molar mass(es) of the system (`[Mw] = kg mol⁻¹`).
"""
get_Mw(eos::Any) = fun_na_error("get_Mw",typeof(eos))

"""
    get_m(eos)

Segment parameter of SAFT-type EOS.
"""
get_m(eos::Any) = fun_na_error("get_m",typeof(eos))

"""
    get_components(eos)

Number of components in the system.
"""
get_components(eos::Any) = fun_na_error("get_components",typeof(eos))