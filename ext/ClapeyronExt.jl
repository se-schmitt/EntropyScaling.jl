module ClapeyronExt

using EntropyScaling, Clapeyron
const ES = EntropyScaling
const CL = Clapeyron

# Bulk properties
ES.pressure(eos::EoSModel, ϱ, T, z=[1.]) = pressure(eos, 1.0./ϱ, T, z)

ES.molar_density(eos::EoSModel, p, T, z=[1.]; phase=:unknown) = molar_density(eos, p, T, z; phase=phase)

ES.entropy_conf(eos::EoSModel, ϱ, T, z=[1.]) = CL.VT_entropy_res(eos, 1.0./ϱ, T, z)

ES.second_virial_coefficient(eos::EoSModel, T, z=[1.]) = second_virial_coefficient(eos, T, z)

# Critical properties
ES.crit_pure(eos::EoSModel) = crit_pure(eos)[1:2]

# Utility functions
ES.split_model(eos::EoSModel) = split_model(eos)

ES.get_Mw(eos::EoSModel) = eos.params.Mw.values .* 1e-3

ES.get_m(eos::EoSModel) = :segment in fieldnames(typeof(eos.params)) ? eos.params.segment.values : ones(length(eos.name))

ES.get_components(eos::EoSModel) = eos.components

end