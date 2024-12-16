#show methods for AbstractEntropyScalingParams
function Base.show(io::IO,params::AbstractEntropyScalingParams)
    print(io,typeof(model).name.name)
    print(io,"(")
    print(io,symbol(params.base.prop))
    print(io,") with fields ")
    print(io,join(fieldnames(typeof(params)),", "))
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
    np = length(model.params)
    for i in 1:length(model.params)
        print(io,name(transport_property(model.params[i])))
        i != np && print(io,", ")
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
    if 
    idx = findfirst(Base.Fix1(transport_compare_type,P),x)
    if transport_compare_type(T,P)
        return quote
            getindex_prop_error(P)
        end
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
            getindex_prop_error(P)
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

#reduced entropy
function reduced_entropy(param::AbstractEntropyScalingParams, s, z=[1.])
    return -s / R / _dot(get_m(param),z)
end

