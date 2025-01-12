# Load Parameters
function load_params(file::String, components::Vector{String}; ref="", ref_id="")
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
        # comps_str = join(components[isnothing.(i_components)],"', '")
        # _file = file[length(DB_PATH)+2:end]
        # msg = "Components ['$comps_str'] not found in database '$_file'!"
        # @info msg
        return missing
    else
        refs_short = unique(String.(data[i_components,j_ref]))
        refs_id = unique(String.(data[i_components,j_ref_id]))
        refs = [Reference(id,short) for (id,short) in zip(refs_id,refs_short)]

        params = [Float64.(data[i_components,j]) for j in j_params]
    
        return params..., refs
    end
end

function load_params(MODEL::Type{<:AbstractTransportPropertyModel}, prop, components; ref="", ref_id="")
    db_path = get_db_path(MODEL, prop)
    return load_params(db_path, lowercase.(components); ref=ref, ref_id=ref_id)
end

function get_db_path(MODEL::Type{<:AbstractTransportPropertyModel}, prop)
    fn = replace(string(MODEL),"Model"=>"")*"_"*replace(name(prop)," "=>"_")*".csv"
    return normpath(DB_PATH, fn)
end
get_db_path(::Type{ChapmanEnskogModel}, prop) = normpath(DB_PATH, "ChapmanEnskog.csv")