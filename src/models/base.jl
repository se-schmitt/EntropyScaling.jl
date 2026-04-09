struct ParamVector{ETA,LAMBDA,D}
    viscosity::ETA
    thermal_conductivity::LAMBDA
    diffusion::D
end
ParamVector(d::Dict) = ParamVector(d[Viscosity()], d[ThermalConductivity()], d[DiffusionCoefficient()])
_get_empty_params_dict() = Dict{AbstractTransportProperty,Any}(prop => missing for prop in [Viscosity(), ThermalConductivity(), DiffusionCoefficient()])
Base.getindex(x::ParamVector, ::Viscosity) = x.viscosity
Base.getindex(x::ParamVector, ::ThermalConductivity) = x.thermal_conductivity
Base.getindex(x::ParamVector, ::AbstractDiffusionCoefficient) = x.diffusion
Base.ismissing(::ParamVector{E,L,D}) where {E,L,D} = _is_Missing(E) && _is_Missing(L) && _is_Missing(D)
_is_Missing(x) = x == Missing
Base.length(::ParamVector{E,L,D}) where {E,L,D} = _is_Missing(E) + _is_Missing(L) + _is_Missing(D)
get_props(params::ParamVector) = begin
    props = AbstractTransportProperty[]
    for field in fieldnames(ParamVector)
        val = getproperty(params, field)
        if !ismissing(val)
            push!(props, transport_property(val))
        end
    end
    return props
end

#show methods for AbstractEntropyScalingParam
function Base.show(io::IO,params::AbstractEntropyScalingParam{P}) where P
    print(io,typeof(params).name.name)
    print(io,"{")
    print(io,symbol(P))
    print(io,"}(;")
    print(io,join(fieldnames(typeof(params)),", "))
    print(io,")")
end

function Base.show(io::IO,::MIME"text/plain",params::AbstractEntropyScalingParam{P}) where P
    print(io,typeof(params).name.name)
    print(io,"{")
    print(io,symbol(P))
    print(io,"} ($(length(params)) components) with fields: ")
    print(io,join(fieldnames(typeof(params)),", "))
end

# show methods for CE model
function Base.show(io::IO, model::ChapmanEnskogModel)
    print(io,typeof(model).name.name)
    print(io,"{")
    print(io,join(model.components,','))
    print(io,"}")
end

function Base.show(io::IO,::MIME"text/plain", model::ChapmanEnskogModel)
    print(io,"ChapmanEnskog{$(join(model.components,','))}")
    print(io,"\n σ: [$(join(round.(model.sigma.values/1e-10,digits=5),", "))] Å")
    print(io,"\n ε: [$(join(round.(model.epsilon.values/kB,digits=3),", "))] K")
    print(io,"\n M: [$(join(round.(model.Mw.values,digits=5),", "))] g/mol")
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
    print(io,join(name.(get_props(model.params)), ", "))
    println(io)
    print(io," Equation of state: ")
    print(io,model.eos)
end

function Base.show(io::IO,model::AbstractEntropyScalingModel)
    print(io,typeof(model).name.name)
    print(io,"(")
    print(io,typeof(model.eos))
    print(io,", ")
    types = map(x -> transport_property(x),model.params)
    print(io,types)
end

#transport_property methods
transport_property(x::AbstractTransportProperty) = x
transport_property(x::BaseParam) = x.prop
transport_property(::AbstractEntropyScalingParam{P}) where P = P()

#Base.getindex methods

function Base.getindex(model::AbstractEntropyScalingModel, prop::P) where P <: AbstractTransportProperty
    return getindex_prop(model.params,prop)
end

function getindex_prop(params::ParamVector, prop::AbstractTransportProperty)
    return params[prop]
end

get_prop_type(m) = typeof(transport_property(m))
get_prop_type(m::Type{BaseParam{P}}) where P = P

function get_prop_type(::Type{T}) where T <: AbstractEntropyScalingParam
    return get_prop_type(fieldtype(T,:base))
end

#valid for any container of AbstractEntropyScalingParam
function getindex_prop(x,prop::P) where P <: AbstractTransportProperty
    idx = findfirst(Base.Fix1(transport_compare_type,P),x)
    if isnothing(idx)
        return getindex_prop_error(prop)
    else
        return x[idx]
    end
end

#when the param is just a AbstractEntropyScalingParam
function getindex_prop(x::AbstractEntropyScalingParam,prop::P) where P <: AbstractTransportProperty
    T = get_prop_type(x)
    if !transport_compare_type(T,P)
        return getindex_prop_error(P)
    else
        return x
    end
end

#=
when the param is a tuple of AbstractEntropyScalingParam.

note to developers,

this function allows access to the properties without allocations, but it is a generated function
so no function inside can be overloaded during a julia session.
=#
@generated function getindex_prop(x::T,prop::P) where {T<:NTuple{<:Any,AbstractEntropyScalingParam},P<:AbstractTransportProperty}
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
    throw(error("cannot find specified property $P"))
end

#get_m

get_m(m::AbstractEntropyScalingParam) = m.m

# entropy scaling variable
function scaling_variable(param::AbstractEntropyScalingParam, s, z = Z1)
    return -s / R
end

function build_model(::Type{MODEL},eos,param_dict::Dict{P}) where {MODEL<:AbstractEntropyScalingModel,P<:AbstractTransportProperty}
    params_vec = Any[]
    PARAM = paramstype(MODEL)
    for (prop, ps) in param_dict
        push!(params_vec, build_param(PARAM, prop, eos, ps))
    end
    params = tuple(params_vec...)
    return MODEL(eos, params)
end

build_model(MODEL,components,params,eos) = MODEL(components,params,eos,String[])
build_model(::Type{MODEL},eos,params::Tuple) where MODEL = build_model(MODEL,get_components(eos),params,_eos_cache(eos))
build_model(::Type{MODEL},eos,params::AbstractVector) where MODEL = build_model(MODEL,eos,tuple(params...))
build_model(::Type{MODEL},eos,params::AbstractEntropyScalingParam) where MODEL = build_model(MODEL,get_components(eos),params,_eos_cache(eos))
function build_param(::Type{PARAM},prop::AbstractTransportProperty,eos,ps) where PARAM <: AbstractEntropyScalingParam
    PARAM(prop,eos,ps...)
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
