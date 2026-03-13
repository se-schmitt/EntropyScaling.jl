using EntropyScaling, Clapeyron, GCIdentifier, ChemicalIdentifiers

println(".")
list = [get_groups_from_name("propanol", gcPCPSAFTGroups)]
println("..")
eos_model = HomogcPCPSAFT(list)
println("...")
model = GCESModel(eos_model, list)
println("....")
η = viscosity(model, 1e5, 200.) #in properties

println("η: ", η)

# using EntropyScaling, Clapeyron
# kB = EntropyScaling.kB

# components = [("ethanol", ["CH3" => 1, "CH2" => 1, "OH" => 1])]

# eos = Clapeyron.HomogcPCPSAFT(components)

# σ = [eos.pcpmodel.params.sigma.values[i,i] for i in 1:Int(sqrt(length(eos.pcpmodel.params.sigma.values)))]
# ϵ = kB .* [eos.pcpmodel.params.epsilon.values[i,i] for i in 1:Int(sqrt(length(eos.pcpmodel.params.epsilon.values)))]
# Mw = eos.pcpmodel.params.Mw.values * 1e-3


# #for i in eachindex(components)
# #    CE_m = EntropyScaling.ChapmanEnskogModel(first.(components)[i], σ[i], ϵ[i], Mw[i])
# #    println(viscosity(CE_m, 300))
# #end    

# #println(eos.params.sigma["CH3"])
# out = EntropyScaling.load_params(EntropyScaling.GCESModel, Viscosity(), components, GC=true)
# A,B,C,D = out
# println(A)
