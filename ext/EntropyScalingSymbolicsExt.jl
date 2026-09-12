module EntropyScalingSymbolicsExt

using EntropyScaling, Symbolics
using Symbolics.SymbolicUtils: BasicSymbolic

const ES = EntropyScaling

@register_symbolic ES.VT_viscosity(model::ES.AbstractEntropyScalingModel, V, T, z::AbstractVector)
@register_symbolic ES.VT_thermal_conductivity(model::ES.AbstractEntropyScalingModel, V, T, z::AbstractVector)
@register_symbolic ES.VT_self_diffusion_coefficient(model::ES.AbstractEntropyScalingModel, V, T, z::AbstractVector)
@register_symbolic ES.VT_self_diffusion_coefficient(model::ES.AbstractEntropyScalingModel, V, T) false
@register_symbolic ES.VT_MS_diffusion_coefficient(model::ES.AbstractEntropyScalingModel, V, T, z::AbstractVector)

for TV in (Num, BasicSymbolic, Real), TT in (Num, BasicSymbolic, Real)
    TV === Real && TT === Real && continue
    @eval ES.VT_thermal_conductivity(model::ES.RefpropRES, V::$TV, T::$TT, z::AbstractVector) =
        Symbolics.wrap(Symbolics.term(ES.VT_thermal_conductivity, model, V, T, z))
end

@register_symbolic ES.viscosity(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)
@register_symbolic ES.viscosity(model::ES.ChapmanEnskogModel, T) false
@register_symbolic ES.thermal_conductivity(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)
@register_symbolic ES.thermal_conductivity(model::ES.ChapmanEnskogModel, T) false
@register_symbolic ES.self_diffusion_coefficient(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)
@register_symbolic ES.self_diffusion_coefficient(model::ES.ChapmanEnskogModel, T) false
@register_symbolic ES.MS_diffusion_coefficient(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)

end
