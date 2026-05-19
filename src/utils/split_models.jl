# Singleton/metadata types are passed through unchanged
CL.is_splittable(::AbstractTransportProperty) = false
CL.is_splittable(::AbstractCollisionIntegralMethod) = false

const SPLIT_TYPES = Union{AbstractEntropyScalingModel,AbstractParam,ChapmanEnskogModel,AbstractParamVector}

function CL._each_split_model(param::T, field, fieldname, group, Ic, Ig) where {T<:SPLIT_TYPES}
    if !CL.is_splittable(field) || fieldname == :sources
        return field
    elseif isnothing(group) || field isa CL.EoSModel
        return CL.each_split_model(field, Ic)
    else
        return CL.each_split_model(field, group, Ic, Ig)
    end
end

CL.default_splitter(model::AbstractEntropyScalingModel) = 1:length(model)

function CL.each_split_model(model::AbstractEntropyScalingModel, I)
    if I isa AbstractVector{Bool}
        return CL.each_split_model(model, findall(I))
    end
    if hasfield(typeof(model), :groups)
        groups = model.groups
        Ig = CL.create_group_splitter(groups, I)
        return CL.each_split_model_struct(model, groups, I, Ig)
    end
    return CL.each_split_model_struct(model, I)
end

CL.each_split_model(model::ChapmanEnskogModel, I) = CL.each_split_model_struct(model, I)
CL.each_split_model(params::AbstractParamVector, I) = CL.each_split_model_struct(params, I)
CL.each_split_model(param::AbstractParam, I) = CL.each_split_model_struct(param, I)

# Vector of property-specific params: split each entry (component-wise)
CL.each_split_model(v::AbstractVector{<:AbstractParam}, I) =
    map(p -> CL.each_split_model(p, I), v)
