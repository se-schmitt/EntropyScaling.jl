second_virial_coefficient_dT(eos::Any, T, z=Z1) = ForwardDiff.derivative(xT -> CL.second_virial_coefficient(eos, xT, z), T)

get_m(eos::Any) = fun_na_error("get_m",typeof(eos))
get_components(eos::Any) = fun_na_error("get_components",typeof(eos))



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
