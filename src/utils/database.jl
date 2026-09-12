function get_db_path(MODEL::Type{<:AbstractTransportPropertyModel}, prop, eos)
    fn = replace(db_model_path(MODEL), "[PROP]" => replace(name(prop)," "=>"_"), "[EOS]" => nameof(typeof(eos)))
    return normpath(DB_PATH, fn)
end

name(model::CL.EoSModel) = split(string(typeof(model)),"{")[1]
name(x) = ""

CL.promote_model(::Type{T}, model::ChapmanEnskogModel) where T = model
CL.promote_model(::Type{T}, prop::AbstractTransportProperty) where T = prop