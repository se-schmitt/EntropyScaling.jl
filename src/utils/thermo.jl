second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> CL.second_virial_coefficient(eos, xT, z), T)

get_components(eos::Any) = fun_na_error("get_components",typeof(eos))
get_components(eos::CL.EoSModel) = eos.components

function CL.init_model(model::AESM, components, userlocations=String[], verbose=false, reference_state=nothing)
    return model
end
function CL.init_model(::Type{𝕄}, components, userlocations=String[], verbose=false, reference_state=nothing) where {𝕄 <: AESM}
    return 𝕄(components)
end
