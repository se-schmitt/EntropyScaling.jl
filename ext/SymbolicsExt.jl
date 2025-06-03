module SymbolicsExt

using EntropyScaling, Symbolics

const ES = EntropyScaling

@register_symbolic ES.ϱT_viscosity(model::ES.AbstractEntropyScalingModel, ϱ, T, z::AbstractVector)
@register_symbolic ES.ϱT_thermal_conductivity(model::ES.AbstractEntropyScalingModel, ϱ, T, z::AbstractVector)
@register_symbolic ES.ϱT_thermal_conductivity(model::ES.RefpropRESModel, ϱ, T, z::AbstractVector)
@register_symbolic ES.ϱT_self_diffusion_coefficient(model::ES.AbstractEntropyScalingModel, ϱ, T, z::AbstractVector)
@register_symbolic ES.ϱT_self_diffusion_coefficient(model::ES.FrameworkModel, ϱ, T, z::AbstractVector)
@register_symbolic ES.ϱT_MS_diffusion_coefficient(model::ES.AbstractEntropyScalingModel, ϱ, T, z::AbstractVector)

@register_symbolic ES.viscosity(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)
@register_symbolic ES.viscosity(model::ES.ChapmanEnskogModel, T)
@register_symbolic ES.thermal_conductivity(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)
@register_symbolic ES.thermal_conductivity(model::ES.ChapmanEnskogModel, T)
@register_symbolic ES.self_diffusion_coefficient(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)
@register_symbolic ES.self_diffusion_coefficient(model::ES.ChapmanEnskogModel, T)
@register_symbolic ES.MS_diffusion_coefficient(model::ES.ChapmanEnskogModel, p, T, z::AbstractVector)

end