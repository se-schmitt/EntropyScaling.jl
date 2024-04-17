module ClapeyronExt

# Extension providing wrappers for functions `call_entropy_scaling` and `fit_entropy_scaling` with the `Clapeyron.jl` (www.github.com/ClapeyronThermo/Clapeyron.jl) package

using EntropyScaling, Clapeyron
const ES = EntropyScaling
const CL = Clapeyron

# Wrapper for the function `fit_entropy_scaling`
```
Test documentation
```
function ES.fit_entropy_scaling(model::EoSModel,
                                T::Vector{Float64}, 
                                ϱ::Vector{Float64}, 
                                Y::Vector{Float64}, 
                                prop::String;
                                i_fit=[0,1,1,1,1],
                                solute::Dict{Symbol,Float64}=Dict{Symbol,Float64}())

    # Creae function handles
    fs = get_Clapeyron_fun(model)

    # Calculate critical temperature and pressure
    (Tc,pc,vc) = crit_pure(model)

    m = model isa SAFTModel ? model.params.segment.values : ones(length(model.name))
    M = model.params.Mw./1e3

    modelDict = Dict(:sfun=>fs.sfun, :Bfun=>fs.Bfun[1], :Tc=>Tc[1], :pc=>pc[1], :M=>M[1], :m_EOS=>m[1])

    return fit_entropy_scaling(modelDict, T, ϱ, Y, prop; i_fit=i_fit, reduced=false, solute=solute)
end


# Wrapper for the function `call_entropy_scaling`
function ES.call_entropy_scaling(model::EoSModel,
                                T::Vector{Float64}, 
                                ϱ::Vector{Float64}, 
                                prop::String; 
                                α::Vector{Vector{Float64}},
                                x::Matrix{Float64}=ones(length(T),1),
                                difcomp::Int64=0)

    # Creae function handles
    fs = get_Clapeyron_fun(model)

    # Calculate critical temperature and pressure
    datc = crit_pure.(split_model(model))
    Tc = [datc[i][1] for i in eachindex(datc)]
    pc = [datc[i][2] for i in eachindex(datc)]

    m = model isa SAFTModel ? model.params.segment.values : ones(length(model.name))
    M = model.params.Mw.values./1e3

    modelDict = Dict(:sfun=>fs.sfun, :Bfun=>fs.Bfun, :Bmixfun=>fs.Bmixfun, :Tc=>Tc, :pc=>pc, :M=>M, :m_EOS=>m, :α=>α)

    return call_entropy_scaling(modelDict, T, ϱ, prop; x=x, reduced=false, difcomp=difcomp)
end

# Function to provide function objectives from Clapeyron
function get_Clapeyron_fun(model)
    Ms = model.params.Mw.values.*1e-3

    # Configurational entropy
    sfun(T,ϱ,x) = [CL.VT_entropy_res(model,x[i,:]'*Ms/ϱ[i],T[i],x[i,:]) for i in eachindex(T)][T isa Vector ? (:) : 1]

    # Second virial coefficient
    Bfun = [T -> second_virial_coefficient.(model_pure, T) for model_pure in split_model(model)]

    # Second virial coefficient of mixture
    Bmixfun(T,x) = [second_virial_coefficient(model, T[i], x[i,:]) for i in eachindex(T)][T isa Vector ? (:) : 1]
    
    return (;sfun=sfun, Bfun=Bfun, Bmixfun=Bmixfun)
end

end