second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> CL.second_virial_coefficient(eos, xT, z), T)

get_components(eos::CL.EoSModel) = eos.components
