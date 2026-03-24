# using EntropyScaling, Clapeyron, GCIdentifier, ChemicalIdentifiers, DelimitedFiles, Plots, JSON
# using Statistics

# # @testset "GCES" begin
# # 	@testset "Pure Viscosity" begin

# function test()
# 	# Test feos

# 	T_k_to_test = [200, 600]
# 	P_MPa_to_test = [1e3, 2e6]

# 	# components_to_test = [("Pentadecane", ["CCCCCCCCCCCCCCC";{"CH3":2;"CH2":13}]), ("Hexane", ["CCCCCC";{"CH3":2;"CH2":4}]), ("Heptane", ["CCCCCCC";{"CH3":2;"CH2":5}]), ("Hexadecane", ["CCCCCCCCCCCCCCCC";{"CH3":2;"CH2":14}]), ("Diallyl ether", ["C=CCOCC=C";{"CH2=":2;"CH=":2;"OCH2":1;"CH2":1}]), ("Monoethanolamine", ["NCCO";{"CH2":2;"OH":1;"NH2":1}]), ("N-Propylamine", ["CCCN";{"CH3":1;"CH2":2;"NH2":1}])]

# 	components_to_test = [
# 		("Pentadecane", [("CCCCCCCCCCCCCCC", ["CH3"=>2, "CH2"=>13])], [1000, 200, 137.2854632919952 * 1e-3]),
# 		("Hexane", [("CCCCCC", ["CH3"=>2, "CH2"=>4])], [2000000, 200, 2.2832750470577547 * 1e-3]),
# 		("Heptane", [("CCCCCCC", ["CH3"=>3, "CH2"=>5])], [2000000, 600, 26.189189122692195 * 1e-6]),
# 	]

# 	for component in components_to_test
# 		list = component[2]
#         println(list)
# 		eos_model = HomogcPCPSAFT(list)
# 		model = GCESModel(eos_model, list)
# 		# @test viscosity(model_ace, P, T)/1e-6 ≈ component[3][i] rtol=1e-4
# 		# println("Test passed: ", viscosity(model, component[3][1], component[3][2])/1e-6 - component[3][3]<1e-4)
# 		println("Test passed: ", viscosity(model, component[3][1], component[3][2]) - component[3][3]<1e-4)
# 		println("difference: ", viscosity(model, component[3][1], component[3][2]) - component[3][3])
# 		println("difference: ",component[3][3])
# 		println("difference: ", viscosity(model, component[3][1], component[3][2]))
# 	end
# end

# test()







@testset "GCES" begin
	@testset "Pure Viscosity" begin
		# Test feos

		T_k_to_test = [200, 600]
		P_MPa_to_test = [1e3, 2e6]
		# components_to_test = [("Pentadecane", ["CCCCCCCCCCCCCCC";{"CH3":2;"CH2":13}]), ("Hexane", ["CCCCCC";{"CH3":2;"CH2":4}]), ("Heptane", ["CCCCCCC";{"CH3":2;"CH2":5}]), ("Hexadecane", ["CCCCCCCCCCCCCCCC";{"CH3":2;"CH2":14}]), ("Diallyl ether", ["C=CCOCC=C";{"CH2=":2;"CH=":2;"OCH2":1;"CH2":1}]), ("Monoethanolamine", ["NCCO";{"CH2":2;"OH":1;"NH2":1}]), ("N-Propylamine", ["CCCN";{"CH3":1;"CH2":2;"NH2":1}])]
		components_to_test = [
			("Pentadecane", [("CCCCCCCCCCCCCCC", ["CH3"=>2, "CH2"=>13])], [1000, 200, 137.2854632919952 * 1e-3]),
			("Hexane", [("CCCCCC", ["CH3"=>2, "CH2"=>4])], [2000000, 200, 2.2832750470577547 * 1e-3]),
			("Heptane", [("CCCCCCC", ["CH3"=>3, "CH2"=>5])], [2000000, 600, 26.189189122692195 * 1e-6]),
		]

		for component in components_to_test
			list = component[2]
			println(list)
			eos_model = HomogcPCPSAFT(list)
			model = GCESModel(eos_model, list)
			@test viscosity(model, component[3][1], component[3][2]) ≈ component[3][3] rtol=1e-4
		end
	end
end
