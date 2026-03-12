using EntropyScaling, Clapeyron, GCIdentifier, ChemicalIdentifiers


list = [get_groups_from_name("propanol", gcPCPSAFTGroups)]

eos_model = HomogcPCPSAFT(list)

model = GCESModel(eos_model, list)

η = viscosity(model, 1e5, 323.) #in properties

println("η: ", η)