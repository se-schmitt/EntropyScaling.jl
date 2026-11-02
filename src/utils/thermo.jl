second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> CL.second_virial_coefficient(eos, xT, z), T)

get_components(eos::CL.EoSModel) = eos.components

# Clapeyron potential models describe one pure fluid with σ = ε = m = 1, no parameters
get_components(eos::CL.PotentialModel) = [string(nameof(typeof(eos)))]
_mw(eos) = CL.mw(eos)
_mw(::CL.PotentialModel) = [1.]
_split(eos, I...) = CL.split_model(eos, I...)
_split(eos::CL.PotentialModel, I...) = [eos]
