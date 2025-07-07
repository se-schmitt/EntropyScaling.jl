module JutulDarcyExt

using EntropyScaling, JutulDarcy

const ES = EntropyScaling
const JD = JutulDarcy

@jutul_secondary function JD.update_viscosity!(mu, es_model::AbstractTransportPropertyModel, 
    model::SimulationModel{D,S}, p, T, flash_results, ix) where {D, S<:JD.CompositionalSystem}

    @inbounds for i in ix
        mu[i] = viscosity(es_model, p, T)
    end

    
    return nothing
end

end