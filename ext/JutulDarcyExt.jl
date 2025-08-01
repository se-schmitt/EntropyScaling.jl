module JutulDarcyExt

using JutulDarcy, Jutul, EntropyScaling

const ES = EntropyScaling
const JD = JutulDarcy

struct JDEntropyScalingModel{M} <: JD.PhaseVariables
    model::M
end
ES.EStoJD(model::ES.AbstractEntropyScalingModel) = JDEntropyScalingModel(model)

@jutul_secondary function update_viscosity!(mu, es::JDEntropyScalingModel{M}, 
    model::JD.SimulationModel, p, T, flash_results, ix) where {M<:ES.AbstractEntropyScalingModel}

    sys = model.system
    l, v = phase_indices(sys)
    @inbounds for i in ix
        mu[l,i] = viscosity(es.model, p, T, flash_results[i])
        mu[v,i] = viscosity(es.model, p, T, flash_results[i])
    end
    
    return nothing
end

end