
@testset "REFPROP Viscosity" begin
    eos_model = IAPWS95()
    #a1,a2,a3,a4,ξ, from https://doi.org/10.1007/s10765-022-03096-9, table 1
    ni = [-0.669333;1.418565;-0.797172;0.160096;;]
    σ = [0.2640e-9]
    ε = [809.100]
    Mw = [18.015e-3]
    model = RefpropRESModel(eos_model,[
        RefpropRESParams(Viscosity(), ni, σ, ε, Mw)
    ])
    #= 
    julia> PropsSI("viscosity","P",1e5,"T",298.15,"water")
    0.0008900226737731678
    =#
    #TODO compare to python
    @test viscosity(model,1e5,298.15) ≈ 0.00089 rtol = 1e-2
end