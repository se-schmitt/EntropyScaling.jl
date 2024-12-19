#show methods for AbstractEntropyScalingParams
function Base.show(io::IO,params::AbstractEntropyScalingParams)
    print(io,typeof(params).name.name)
    print(io,"{")
    print(io,symbol(params.base.prop))
    print(io,"}(;")
    print(io,join(fieldnames(typeof(params)),", "))
    print(io,")")
end

function Base.show(io::IO,::MIME"text/plain",params::AbstractEntropyScalingParams)
    print(io,typeof(params).name.name)
    print(io,"{")
    print(io,symbol(params.base.prop))
    print(io,"} ($(length(params.base.Mw)) components) with fields: ")
    print(io,join(fieldnames(typeof(params)),", "))
end

#show methods for CE model
function Base.show(io::IO, model::ChapmanEnskogModel)
    print(io,typeof(params).name.name)
    print(io,"{")
    print(io,join(model.components,','))
    print(io,"}")
end

function Base.show(io::IO,::MIME"text/plain", model::ChapmanEnskogModel)
    print(io,"ChapmanEnskogModel{$(join(model.components,','))}")
    print(io,"\n σ: [$(join(round.(model.σ/1e-10,digits=5),", "))] Å")
    print(io,"\n ε: [$(join(round.(model.ε/kB,digits=3),", "))] K")
    print(io,"\n M: [$(join(round.(model.Mw,digits=5),", "))] kg/m³")
    print(io,"\n Collision integral: $(typeof(model.collision).name.name)")
end

#show methods for AbstractEntropyScalingModel
function Base.show(io::IO,::MIME"text/plain",model::AbstractEntropyScalingModel)
    n = length(model.components)
    print(io,typeof(model).name.name)
    print(io," with ",n," component")
    n != 1 && print(io,"s")
    println(io,":")
    Base.print_matrix(IOContext(io, :compact => true),model.components)
    println(io)
    print(io," Available properties: ")
    if model.params isa AbstractEntropyScalingParams
        print(io,name(transport_property(model.params)))
    else
        np = length(model.params)
        for i in 1:length(model.params)
            print(io,name(transport_property(model.params[i])))
            i != np && print(io,", ")
        end
    end
    println(io)
    print(io," Equation of state: ")
    print(io,model.eos)
end

function Base.show(io::IO,model::AbstractEntropyScalingModel)
    print(io,typeof(model).name.name)
    print(io,"(")
    print(io,typeof(model.eos))
    print(io,", ")
    types = map(x -> x.base.prop,model.params)
    print(io,types)
end

#transport_property methods
transport_property(x::AbstractTransportProperty) = x
transport_property(x::BaseParam) = x.prop
transport_property(x::AbstractEntropyScalingParams) = transport_property(x.base)

#Base.getindex methods

function Base.getindex(model::AbstractEntropyScalingModel, prop::P) where P <: AbstractTransportProperty
    return getindex_prop(model.params,prop)
end

get_prop_type(m) = typeof(transport_property(m))
get_prop_type(m::Type{BaseParam{P}}) where P = P

function get_prop_type(::Type{T}) where T <: AbstractEntropyScalingParams
    return get_prop_type(fieldtype(T,:base))
end

#valid for any container of AbstractEntropyScalingParams
function getindex_prop(x,prop::P) where P <: AbstractTransportProperty
    idx = findfirst(Base.Fix1(transport_compare_type,P),x)
    if isnothing(idx)
        return quote
            getindex_prop_error(P)
        end
    else
        return x[idx]
    end
end

#when the param is just a AbstractEntropyScalingParam
function getindex_prop(x::AbstractEntropyScalingParams,prop::P) where P <: AbstractTransportProperty
    T = get_prop_type(x)
    if !transport_compare_type(T,P)
        return getindex_prop_error(P)
    else
        return x
    end
end

#=
when the param is a tuple of AbstractEntropyScalingParams.

note to developers,

this function allows access to the properties without allocations, but it is a generated function
so no function inside can be overloaded during a julia session.
=#
@generated function getindex_prop(x::T,prop::P) where {T<:NTuple{<:Any,AbstractEntropyScalingParams},P<:AbstractTransportProperty}
    idx = findfirst(xi -> transport_compare_type(get_prop_type(xi),P),fieldtypes(T))
    if isnothing(idx)
        return quote
            getindex_prop_error(prop)
        end
    else
        f = fieldtypes(T)[idx]
        return :(x[$idx]::$(f))
    end
end

function getindex_prop_error(P)
    throw(error("cannot found specified property $P"))
end

#get_m

get_m(m::AbstractEntropyScalingParams) = m.m

#caching
_eos_cache(eos) = eos

# entropy scaling variable
function scaling_variable(param::AbstractEntropyScalingParams, s, z = Z1)
    return -s / R
end

function build_model(::Type{MODEL},eos,param_dict::Dict{P}) where {MODEL<:AbstractEntropyScalingModel,P<:AbstractTransportProperty}
    params_vec = Any[]
    PARAM = paramstype(MODEL)
    for (prop, α) in param_dict
        T = eltype(α)
        if α isa Vector && length(eos) == 1
            αx = convert(Array{T,2},reshape(α,length(α),1))
        else
            αx = convert(Array{T,2},α)
        end
        push!(params_vec, build_param(PARAM, prop, eos, αx))
    end
    params = tuple(params_vec...)
    return MODEL(eos, params)
end

build_model(MODEL,components,params,eos) = MODEL(components,params,eos)
build_model(::Type{MODEL},eos,params::Tuple) where MODEL = build_model(MODEL,get_components(eos),params,_eos_cache(eos))
build_model(::Type{MODEL},eos,params::AbstractVector) where MODEL = build_model(MODEL,eos,tuple(params...))
build_model(::Type{MODEL},eos,params::AbstractEntropyScalingParams) where MODEL = build_model(MODEL,get_components(eos),params,_eos_cache(eos))
function build_param(::Type{PARAM},prop::AbstractTransportProperty,eos,αx;kwargs...) where PARAM <: AbstractEntropyScalingParams
    PARAM(prop,eos,αx,kwargs...)
end

function paramstype end
#a macro that automates some definitions of model methods.
#TODO: handle fitting for arbitrary methods.
macro modelmethods(model,param)
    return quote
        EntropyScaling.paramstype(::Type{<:$model}) = $param
        $model(eos,params::Tuple) = EntropyScaling.build_model($model,eos,params)
        $model(eos,params::AbstractVector) = EntropyScaling.build_model($model,eos,params)
        $model(eos,params::Dict{P}) where P<:EntropyScaling.AbstractTransportProperty = EntropyScaling.build_model($model,eos,params)
        $model(eos,params::$param) = EntropyScaling.build_model($model,eos,params)
    end |> esc
end
