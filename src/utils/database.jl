# Load Parameters
function load_params(file::String, components::Vector{<:AbstractString}; ref="", ref_id="")
    data, header = readdlm(file, ','; header=true)
    j_subs = findfirst(header[:] .== "substance") 
    j_ref = findfirst(header[:] .== "ref")
    j_ref_id = findfirst(header[:] .== "ref_id")
    @assert !any(isnothing.([j_subs,j_ref,j_ref_id])) "Wrong format of header in file $file." 
    N_cols = length(header)
    j_params = findall((!).(in).(1:N_cols,Ref(vcat(j_subs,j_ref,j_ref_id))))
    
    what_refs = contains.(data[:,j_ref][:],ref) .&& contains.(data[:,j_ref_id][:],ref_id)
    subs = split.(data[:,j_subs][:],"|")
    i_components = [findfirst(in.(c,subs) .&& what_refs) for c in components]
    if any(isnothing.(i_components))
        return missing
    else
        refs_short = unique(String.(data[i_components,j_ref]))
        refs_id = unique(String.(data[i_components,j_ref_id]))
        refs = [Reference(id,short) for (id,short) in zip(refs_id,refs_short)]

        params = [Float64.(data[i_components,j]) for j in j_params]
    
        return params..., refs
    end
end

function load_params_GC(file::String, components::Vector; ref="", ref_id="")
    data, header = readdlm(file, ','; header=true)
    j_groups = findfirst(header[:] .== "groups") #Spaltenindex für groups
    j_ref = findfirst(header[:] .== "ref")
    j_ref_id = findfirst(header[:] .== "ref_id")
    @assert !any(isnothing.([j_groups,j_ref,j_ref_id])) "Wrong format of header in file $file." 
    N_cols = length(header)
    j_params = findall((!).(in).(1:N_cols,Ref(vcat(j_groups,j_ref,j_ref_id)))) #Spaltenidizes für Parameter

    groups = reduce(vcat, last.(components))
    A_a = Dict()
    B_a = Dict()
    C_a = Dict()
    D_a = Dict()
    for group in groups
        if !haskey(A_a, group[1])
            i_component = findfirst(data[:, j_groups] .== group[1])
            A_a[group[1]] = data[i_component, j_params[1]]
            B_a[group[1]] = data[i_component, j_params[2]]
            C_a[group[1]] = data[i_component, j_params[3]]
            D_a[group[1]] = data[i_component, j_params[4]]
        end
    end
    return [A_a, B_a, C_a, D_a]
end

function load_params(MODEL::Type{<:AbstractTransportPropertyModel}, prop, components; ref="", ref_id="", GC=false)
    db_path = get_db_path(MODEL, prop)
    if !GC
        return load_params(db_path, lowercase.(components); ref=ref, ref_id=ref_id)
    else
        return load_params_GC(db_path, components; ref=ref, ref_id=ref_id)
    end
end

function load_refprop_names(_components)
    components = lowercase.(_components)
    refprop_names = fill("",size(components))
    for prop in [Viscosity(), ThermalConductivity()]
        data, header = readdlm(get_db_path(RefpropRESModel, prop), ','; header=true)
        j_subs = findfirst(header[:] .== "substance")
        subs = split.(data[:,j_subs][:],"|")
        for (i,c) in enumerate(components)
            idx = findfirst(in.(c,subs))
            if !isnothing(idx) && isempty(refprop_names[i])
                refprop_names[i] = String(subs[idx][end])
            end
        end
        !any(isempty.(refprop_names)) && break
    end

    @assert !any(isempty.(refprop_names)) "Model(s) for $(_components[isnothing.(i_components)]) not available!"
    
    return refprop_names
end

function get_db_path(MODEL::Type{<:AbstractTransportPropertyModel}, prop)
    MODEL_str = replace(string(MODEL),"Model"=>"", "EntropyScaling."=>"")
    fn = MODEL_str*"_"*replace(name(prop)," "=>"_")*".csv"
    return normpath(DB_PATH, fn)
end
get_db_path(::Type{<:AbstractChapmanEnskogModel}, prop) = normpath(DB_PATH, "ChapmanEnskog.csv")