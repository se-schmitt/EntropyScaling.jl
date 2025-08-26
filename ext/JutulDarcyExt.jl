module JutulDarcyExt

using JutulDarcy, Jutul, EntropyScaling

const ES = EntropyScaling
const JD = JutulDarcy
const FD = Jutul.ForwardDiff

struct JDEntropyScalingModel{M} <: JD.PhaseVariables
    model::M
end
ES.EStoJD(model::ES.AbstractEntropyScalingModel) = JDEntropyScalingModel(model)

@jutul_secondary function update_viscosity!(mu, es::JDEntropyScalingModel{M}, 
    model::JD.SimulationModel, Pressure, Temperature, PhaseMassDensities, FlashResults, ix) where {M<:ES.AbstractEntropyScalingModel}

    sys = model.system
    l, v = JD.phase_indices(sys)

    @inbounds for i in ix
        Mwl = es.model.params[1].base.Mw' * FlashResults[i].liquid.mole_fractions
        Mwv = es.model.params[1].base.Mw' * FlashResults[i].vapor.mole_fractions
        ϱl = PhaseMassDensities[l,i] / Mwl
        ϱv = PhaseMassDensities[v,i] / Mwv

        mu[l,i] = ES.ϱT_viscosity(es.model, FD.value(ϱl), FD.value(Temperature[i]), FD.value.(FlashResults[i].liquid.mole_fractions))
        mu[v,i] = ES.ϱT_viscosity(es.model, FD.value(ϱv), FD.value(Temperature[i]), FD.value.(FlashResults[i].vapor.mole_fractions))
    end
    
    return nothing
end

end