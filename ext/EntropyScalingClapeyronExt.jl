module EntropyScalingClapeyronExt

using EntropyScaling, Clapeyron
const ES = EntropyScaling
const CL = Clapeyron
const Downloads = CL.Downloads
const SA1 = CL.SA[1.0]
const url_refprop = "https://raw.githubusercontent.com/usnistgov/fastchebpure/50af5c154a113ac27a2c0a1c3538bc4f43a73a66/teqp_REFPROP10/dev/fluids/"

# Bulk properties
ES.pressure(eos::EoSModel, ϱ, T, z = SA1) = pressure(eos, 1.0./ϱ, T, z)

ES.molar_density(eos::EoSModel, p, T, z = SA1; phase=:unknown, ϱ0=nothing) = begin 
    molar_density(eos, p, T, z; phase, vol0=isnothing(ϱ0) ? nothing : 1/ϱ0)
end
ES.molar_density(eos::CL.EoSVectorParam, p, T, z = SA1; phase=:unknown, ϱ0=nothing) = begin
    ES.molar_density(eos.model, p, T, z; phase, ϱ0)
end

#note: On Clapeyron, each EoSModel can set its own gas constant. here we divide by that and multiply 
#by the constant used by EntropyScaling, so calculations are consistent between models
ES.entropy_conf(eos::EoSModel, ϱ, T, z = SA1) = CL.VT_entropy_res(eos, 1.0/ϱ, T, z)*ES.R/CL.Rgas(eos)
ES.entropy_conf(eos::CL.EoSVectorParam, ϱ, T, z = SA1) = ES.entropy_conf(eos.model, ϱ, T, z)

ES.second_virial_coefficient(eos::EoSModel, T, z = SA1) = second_virial_coefficient(eos, T, z)
ES.second_virial_coefficient(eos::CL.EoSVectorParam, T, z = SA1) = ES.second_virial_coefficient(eos.model, T, z)

ES.isobaric_heat_capacity(eos::EoSModel, ϱ, T, z = SA1) = CL.VT_isobaric_heat_capacity(eos, 1.0/ϱ, T, z)
ES.isobaric_heat_capacity(eos::CL.EoSVectorParam, ϱ, T, z = SA1) = ES.isobaric_heat_capacity(eos.model, ϱ, T, z)

ES.isochoric_heat_capacity(eos::EoSModel, ϱ, T, z = SA1) = CL.VT_isochoric_heat_capacity(eos, 1.0/ϱ, T, z)
ES.isochoric_heat_capacity(eos::CL.EoSVectorParam, ϱ, T, z = SA1) = ES.isochoric_heat_capacity(eos.model, ϱ, T, z)

# Critical properties
ES.crit_pure(eos::EoSModel) = begin
    Tc,pc,vc = crit_pure(eos)
    return Tc,pc,1/vc
end
ES.crit_mix(eos::EoSModel, z) = begin
    Tc,pc,vc = crit_mix(eos, z)
    return Tc,pc,1/vc
end

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

# Model specific wrapper 
ES.RefpropRESModel(eos::CL.MultiFluid; ηref=nothing) = RefpropRESModel(eos, eos.components; ηref)
ES.RefpropRESModel(comps::AbstractString; kwargs...) = RefpropRESModel([comps]; kwargs...)
function ES.RefpropRESModel(comps::Vector{<:AbstractString}; ηref=nothing, kwargs_CL...)

    kw_CL = NamedTuple(kwargs_CL)
    eos = nothing
    if !(:estimate_mixing in keys(kw_CL))
        kw_CL = merge(kw_CL, (;estimate_mixing=:lb))
    end

    try 
        eos = MultiFluid(comps; kw_CL...)
    catch e1 
        if (e1 isa ErrorException && any(startswith.(Ref(e1.msg),["Coolprop: key","cannot "]))) || e1 isa MissingException
            _comps = copy(comps)
            try 
                names = uppercase.(ES.load_refprop_names(comps))
                _comps .= Downloads.download.(url_refprop .* names .* ".json") .|> read .|> String
            catch e2
                if e2 isa Downloads.RequestError
                    i_err = findfirst(comps .== _comps)
                    throw(ErrorException("Refprop parameters for '$(comps[i_err])' could not be loaded."))
                else 
                    rethrow(e2)
                end
            end
            mixing_userlocations = :mixing_userlocations in keys(kw_CL) ? kw_CL.mixing_userlocations : String[]
            mixing = CL.init_model(AsymmetricMixing,comps,mixing_userlocations,false)
            eos = MultiFluid(_comps; mixing, coolprop_userlocations=false, kw_CL...)
            eos.components .= comps
        else
            rethrow(e1)
        end
    end 

    return RefpropRESModel(eos, comps; ηref)
end

end #module
