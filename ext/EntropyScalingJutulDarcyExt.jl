module EntropyScalingJutulDarcyExt

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
        Mwl = es.model.params[1].base.Mw.values' * FlashResults[i].liquid.mole_fractions
        Mwv = es.model.params[1].base.Mw.values' * FlashResults[i].vapor.mole_fractions
        Vl = Mwl / FD.value(PhaseMassDensities[l,i])
        Vv = Mwv / FD.value(PhaseMassDensities[v,i])

        mu[l,i] = ES.VT_viscosity(es.model, Vl, FD.value(Temperature[i]), FD.value.(FlashResults[i].liquid.mole_fractions))
        mu[v,i] = ES.VT_viscosity(es.model, Vv, FD.value(Temperature[i]), FD.value.(FlashResults[i].vapor.mole_fractions))
    end
    
    return nothing
end

end