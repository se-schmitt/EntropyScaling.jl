module ClapeyronExt

using EntropyScaling, Clapeyron
const ES = EntropyScaling
const CL = Clapeyron
const SA1 = CL.SA[1.0]
# Bulk properties
ES.pressure(eos::EoSModel, ϱ, T, z = SA1) = pressure(eos, 1.0./ϱ, T, z)

ES.molar_density(eos::EoSModel, p, T, z = SA1; phase=:unknown) = molar_density(eos, p, T, z; phase)
ES.molar_density(eos::CL.EoSVectorParam, p, T, z = SA1; phase=:unknown) = ES.molar_density(eos.model, p, T, z; phase)

#note: On Clapeyron, each EoSModel can set its own gas constant. here we divide by that and multiply 
#by the constant used by EntropyScaling, so calculations are consistent between models
ES.entropy_conf(eos::EoSModel, ϱ, T, z = SA1) = CL.VT_entropy_res(eos, 1.0./ϱ, T, z)*ES.R/CL.Rgas(eos)
ES.entropy_conf(eos::CL.EoSVectorParam, ϱ, T, z = SA1) = ES.entropy_conf(eos.model, ϱ, T, z)

ES.second_virial_coefficient(eos::EoSModel, T, z = SA1) = second_virial_coefficient(eos, T, z)
ES.second_virial_coefficient(eos::CL.EoSVectorParam, T, z = SA1) = ES.second_virial_coefficient(eos.model, T, z)

# Critical properties
ES.crit_pure(eos::EoSModel) = crit_pure(eos)[1:2]

# Utility functions
ES.split_model(eos::EoSModel) = split_model(eos)
ES.split_model(eos::SingleFluid) = [eos]
ES.split_model(eos::MultiFluid) = eos.pures
ES.split_model(eos::CL.EoSVectorParam) = eos.pure
ES.get_Mw(eos::EoSModel) = Clapeyron.mw(eos) .* 1e-3
ES.get_Mw(eos::CL.EoSVectorParam) = ES.get_Mw(eos.model)

ES.get_m(eos::EoSModel) = :segment in fieldnames(typeof(eos.params)) ? eos.params.segment.values : ones(length(eos.name))
ES.get_m(eos::CL.EoSVectorParam) = ES.get_m(eos.model)
ES.get_m(eos::SingleFluid) = SA1
ES.get_m(eos::CL.SAFTgammaMieModel) = ES.get_m(eos.vr_model)
ES.get_components(eos::EoSModel) = eos.components

ES._eos_cache(eos::EoSModel) = CL.EoSVectorParam(eos)
ES._eos_cache(eos::CL.EoSVectorParam) = eos
ES._eos_cache(eos::MultiFluid) = eos

end #module
