module MicthermExt

# Extension providing wrappers for functions `call_entropy_scaling` and `fit_entropy_scaling` to be coupled with the `MicTherm` (LINK) package

using EntropyScaling, MATLAB
const ES = EntropyScaling

# MicTherm parameter struct
struct MicThermParam <: ES.MicThermParamType
    name::Vector{String}
    eos::String
    polar::Union{String,Missing}
    assoc::Union{String,Missing}
    unit::String
    M::Union{Vector{Float64},Missing}
    m::Union{Vector{Float64},Missing}
    σ::Union{Vector{Float64},Missing}
    ε::Union{Vector{Float64},Missing}
    ξ::Union{Vector{Float64},Missing}
    Tcm::Union{Vector{Float64},Missing}
    pcm::Union{Vector{Float64},Missing}
    λr::Union{Vector{Float64},Missing}
    λa::Union{Vector{Float64},Missing}
    TStar::Union{Vector{Float64},Missing}
    vStar::Union{Vector{Float64},Missing}
    μ::Union{Vector{Float64},Missing}
    N_μ::Union{Vector{Float64},Missing}
    Q::Union{Vector{Float64},Missing}
    N_Q::Union{Vector{Float64},Missing}
    type_AB::Union{Vector{String},Missing}
    κ_AB::Union{Vector{Float64},Missing}
    ε_AB::Union{Vector{Float64},Missing}
    α::Union{Vector{Float64},Vector{Vector{Float64}},Missing}
    source::Union{Vector{String},Missing}
    path::String
end

"""
`MicThermParam(ParamDict::Dict{Symbol,Any})`

Constructor for MicThermParam struct

- Mandatory input parameters: `:name`, `:eos`, `:path`, `:unit`
- Optional input parameters: `α[_η/_λ/_D]`, `:polar`, `:assoc`, `:M`, `:m`, `:σ`, `:ε`, `:ξ`, `:Tcm`, `:pcm`, `:λr`, `:λa`, `:TStar`, `:vStar`, `:μ`, `:N_μ`, `:Q`, `:N_Q`, `:type_AB`, `:κ_AB`, `:ε_AB`, `:source`

Examplet for constructing a MicThermParam struct:
```julia-repl
d = Dict(:name=>["hexane", "heptane"], :eos=>"PC_SAFT", :path=>"/path/to/MicTherm", :unit=>"SI", :M=>[86.18, 100.21], :m=>[3.0576,3.4831], :σ=>[3.7983,3.8049], :ε=>[236.77,238.4], :source=>["10.1021/ie0003887"])
model = MicThermParam(d)
```
"""
function ES.MicThermParam(ParamDict_ori::Dict{Symbol,Any}; prop=nothing)
    ParamDict = deepcopy(ParamDict_ori)

    # Check Dict
    mandatory_input = [:name, :eos, :path, :unit]
    for f in mandatory_input
        if !(f in keys(ParamDict))  error("Mandatory input parameter $f is missing.") end
    end

    # Set α
    if !(:α in keys(ParamDict))
        sym = prop == "vis" ? :α_η : prop == "tcn" ? :α_λ : prop == "dif" ? :α_D : nothing
        if !isnothing(prop) && sym in keys(ParamDict)
            ParamDict[:α] = ParamDict[sym]
        end
    end

    fields = fieldnames(MicThermParam)
    return MicThermParam([f in keys(ParamDict) ? ParamDict[f] : missing for f in fields]...)
end

# Wrapper for the function `fit_entropy_scaling` to be used with MicTherm
function ES.fit_entropy_scaling(model::ES.MicThermParamType,
                                T::Vector{Float64}, 
                                ϱ::Vector{Float64}, 
                                Y::Vector{Float64}, 
                                prop::String; 
                                i_fit=[0,1,1,1,1],
                                solute::Dict{Symbol,Float64}=Dict{Symbol,Float64}())
    
    # Create function handles
    fs = get_MicTherm_fun(model)

    # Calculate critical temperature and pressure
    (Tc, pc, ϱc) = calc_crit_MicTherm(model)
    
    m = ismissing(model.m) ? ones(length(model.name)) : model.m
    M = model.unit == "reduced" ? model.M : model.M./1e3

    modelDict = Dict(:sfun=>fs.sfun, :Bfun=>fs.Bfun[1], :dBdTfun=>fs.dBdTfun[1], :Tc=>Tc[1], :pc=>pc[1], :M=>M[1], :m_EOS=>m[1])
    if model.unit == "reduced"
        modelDict[:σ] = model.σ[1]
        modelDict[:ε] = model.ε[1]
    end

    return fit_entropy_scaling(modelDict, T, ϱ, Y, prop; i_fit=i_fit, reduced=model.unit=="reduced", solute=solute)
end

# Wrapper for the function `call_entropy_scaling` to be used with MicTherm
function ES.call_entropy_scaling(model::ES.MicThermParamType,
                                T::Vector{Float64}, 
                                ϱ::Vector{Float64},
                                prop::String; 
                                x::Matrix{Float64}=ones(length(T),1),
                                difcomp::Int64=0)
    
    # Create function handles
    fs = get_MicTherm_fun(model)

    # Calculate critical temperature and pressure
    (Tc, pc, ϱc) = calc_crit_MicTherm(model)

    m = ismissing(model.m) ? ones(length(model.name)) : model.m
    M = model.unit == "reduced" ? model.M : model.M./1e3

    modelDict = Dict(:sfun=>fs.sfun, :Bfun=>fs.Bfun, :dBdTfun=>fs.dBdTfun, :Bmixfun=>fs.Bmixfun, :dBdTmixfun=>fs.dBdTmixfun, :Tc=>Tc, :pc=>pc, :M=>M, :m_EOS=>m, :α=>model.α)
    if model.unit == "reduced"
        modelDict[:σ] = model.σ
        modelDict[:ε] = model.ε
        if !ismissing(model.ξ)
            modelDict[:ξ] = model.ξ[1]
        end
    end

    return call_entropy_scaling(modelDict, T, ϱ, prop; x=x, reduced=model.unit=="reduced", difcomp=difcomp)
end

# Function to create argumentsstring for Initialization
function get_initialization_string(p::MicThermParam; comp=0)
    if comp == 0
        N_comp = length(p.name)
    else
        N_comp = 1
    end
    
    out = "'uninitialized', [], [], [], [], 'N_components = $(N_comp)', 'units = $(p.unit)', 'EOS = $(p.eos)'"

    # Polar and association
    for s in [:polar, :assoc]
        if !ismissing(getfield(p,s))
            out *= ", '$(string(s)) = $(getfield(p,s))'"
        end
    end

    # Component parameters
    comps = comp == 0 ? (1:N_comp) : comp
    for i in comps
        if comp == 0 j = i else j = 1 end
        out *= ", 'Substance_ID$j = 0'$(set_param(p,:name,"PotModel_$j",i))"
        
        # EOS all
        out *= set_param(p,:M,"molar_mass_$j",i)

        # EOS specific
        if lowercase(p.eos) in ["pc_saft","kolafanezbeda","saft_vr_mie","backone","stephan","soft_saft","pets"]
            out *= string(set_param.(Ref(p),[:σ,:ε],["sigma_$j","epsilon_$j"],i)...)
        elseif lowercase(p.eos) in ["scpa","cpa"]
            out *= string(set_param.(Ref(p),[:Tcm,:pcm,:m],["Tcm_$j","pcm_$j","m_VDW_$j"],i)...)
        elseif lowercase(p.eos) in ["pactplusb"]
            out *= string(set_param.(Ref(p),[:vStar,:TStar,:vStar,:TStar],["vStar_$j","TStar_$j","sigma_$j","epsilon_$j"],i)...)
        else
            error("Either (σ, ϵ, m) or (Tcm, pcm) or (vTStar, vStar) should be defined!")
        end

        if lowercase(p.eos) in ["pc_saft","saft_vr_mie","backone","soft_saft"]
            out *= set_param.(Ref(p),:m,"chainlength_$j",i)
        end
        
        if lowercase(p.eos) in ["saft_vr_mie"]
            out *= string(set_param.(Ref(p),[:λr,:λa],["lambda_r_$j","lambda_a_$j"],i)...)
        end

        if lowercase(p.eos) in ["backone"]
            out *= set_param(p,:m,["alpha_$j"],i)
        end

        # Polar contribution
        if !ismissing(p.N_μ)
            out *= string(set_param.(Ref(p),[:μ,:N_μ],["dipolemoment_$j","N_dipolemoment_$j"],i)...)
        end
        if !ismissing(p.N_Q)
            out *= string(set_param.(Ref(p),[:Q,:N_Q],["quadrupolemoment_$j","N_quadrupolemoment_$j"],i)...)
        end

        # association
        if !ismissing(p.type_AB)
            out *= string(set_param.(Ref(p),[:type_AB,:κ_AB,:ε_AB],["assoc_type_$j","assoc_kappa_$j","assoc_epsilon_$j"],i)...)
        end
    end

    # Mixing parameter
    if N_comp == 2 && !ismissing(p.ξ)
        out *= set_param(p,:ξ,"Xsi_12",1)
    end

    return out
end

function set_param(p,field,str,i)
    if ismissing(getfield(p,field))
        error("Field $(field) is missing!")
    else
        return ", '$(string(str)) = $(getfield(p,field)[i])'"
    end
end

# Function to provide function objectives from MicTherm
function get_MicTherm_fun(model)
    # Set MicTherm path
    eval_string("addpath(genpath('$(model.path)'))")
    mat"warning off"

    # Set conversion factors
    if model.unit == "reduced"
        Bconv = 1
        ϱconv = 1
        pconv = 1
        h = 1e-3
    elseif model.unit == "SI"
        Bconv = 1e-3
        ϱconv = model.M
        pconv= 1e6
        h = 5
    end

    # Configurational entropy
    sfun(T,ϱ,x) = 
    (   index = T isa Number ? 1 : 1:length(T);
        Mmix = x * model.M;
        eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = sRes');");
        mat"[$names , ~ , $values] = IO_API( 'initialized', $(ϱ./Mmix), $T, [], $x );";                    
        return values[index,findfirst(names[:] .== "sRes")] )

    # Second virial coefficient
    Bfun = [T ->
        (   index = T isa Number ? 1 : 1:length(T);
            ϱ = T isa Number ? 1.0 : ones(length(T));
            x = T isa Number ? 1.0 : ones(length(T),1);
            eval_string("IO_API( $(get_initialization_string(model; comp=i)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = B');");
            mat"[$names , ~ , $values] = IO_API( 'initialized', $ϱ, $T, [], $x );";                    
            return values[index,findfirst(names[:] .== "B")] .* Bconv 
        )
        for i in eachindex(model.name) ]

    # Second virial coefficient of mixture
    Bmixfun(T,x) = 
    (   index = T isa Number ? 1 : 1:length(T);
        ϱ = T isa Number ? 1.0 : ones(length(T));
        eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = B');");
        mat"[$names , ~ , $values] = IO_API( 'initialized', $ϱ, $T, [], $x );";                    
        return values[index,findfirst(names[:] .== "B")] .* Bconv 
    )
    
    # Derivative of second virial coefficient
    dBdTfun = [T -> 
        (   index = T isa Number ? 1 : 1:length(T);
            ϱ = T isa Number ? 1.0 : ones(length(T));
            x = T isa Number ? 1.0 : ones(length(T),1);
            eval_string("IO_API( $(get_initialization_string(model; comp=i)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = B');");
            mat"[$names , ~ , $values_h] = IO_API( 'initialized', $ϱ, $(T.+h/2), [], $x );";
            mat"[$names , ~ , $values_l] = IO_API( 'initialized', $ϱ, $(T.-h/2), [], $x );";
            B_h = values_h[index,findfirst(names[:] .== "B")];
            B_l = values_l[index,findfirst(names[:] .== "B")];
            return (B_h .- B_l) ./ h .* Bconv
        )
        for i in eachindex(model.name) ]
    
    # Derivative of second virial coefficient of mixture
    dBdTmixfun(T,x) = 
    (   index = T isa Number ? 1 : 1:length(T);
        ϱ = T isa Number ? 1.0 : ones(length(T));
        eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = B');");
        mat"[$names , ~ , $values_h] = IO_API( 'initialized', $ϱ, $(T.+h/2), [], $x );";
        mat"[$names , ~ , $values_l] = IO_API( 'initialized', $ϱ, $(T.-h/2), [], $x );";
        B_h = values_h[index,findfirst(names[:] .== "B")];
        B_l = values_l[index,findfirst(names[:] .== "B")];
        return (B_h .- B_l) ./ h .* Bconv
    )

    # Pressure
    pfun(T,ϱ; x=ones(length(T),1)) = 
    (   index = T isa Number ? 1 : 1:length(T);
        Mmix = x * model.M;
        eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = p');");
        mat"[$names , ~ , $values] = IO_API( 'initialized', $(ϱ./Mmix), $T, [], $x );";                    
        return values[index,findfirst(names[:] .== "p")].*pconv 
    )

    # Density
    ϱfun(T,p; x=ones(length(T),1), states=repeat(["L"],length(T)), ϱmax::Number=Inf) = 
    (   index = T isa Number ? 1 : 1:length(T);
        Mmix = x * model.M;
        eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = rho');");
        mat"[$names , ~ , $values] = IO_API( 'initialized', [], $T, $(p./pconv), $x );";
        ϱ_all = values[:,findfirst(names[:] .== "rho")];
        nr = values[:,findfirst(names[:] .== "Nr")];
        ϱ = NaN*ones(length(T));
        for i in eachindex(T);
            what_i = nr .== i .&& ϱ_all .<= ϱmax/Mmix[i];
            ϱ[i] = states[i] == "L" ? maximum(ϱ_all[what_i]) : minimum(ϱ_all[what_i]);
        end;
        return ϱ[index].*Mmix[index]
    )

    funs = (;sfun=sfun, Bfun=Bfun, dBdTfun=dBdTfun, Bmixfun=Bmixfun, dBdTmixfun=dBdTmixfun, pfun=pfun, ϱfun=ϱfun)
    
    return funs
end 

# Function to calculate critical point by MicTherm
function calc_crit_MicTherm(model)
    # Set conversion factors
    if model.unit == "reduced"
        pconv = 1
        ΔT = 0.1
    elseif model.unit == "SI"
        pconv = 1e6
        ΔT = 100
    end
    ϱconv = model.M

    Tc = Float64[]
    pc = Float64[]
    ϱc = Float64[]
    for i in 1:length(model.name)
        eval_string("IO_API( $(get_initialization_string(model; comp=i)), 'calculationmode = VLE', 'dT = $(ΔT)' );")
        mat"[$names_VLE , ~ , $values_VLE] = IO_API( 'initialized', [], [], [], [] );"
        mat"close all;"
        push!(Tc,values_VLE[1,findfirst(names_VLE[:] .== "T_VLE")])          # K or reduced
        push!(pc,values_VLE[1,findfirst(names_VLE[:] .== "p_VLE")]*pconv)    # Pa or reduced
        push!(ϱc,values_VLE[1,findfirst(names_VLE[:] .== "rhoL")]*ϱconv[i])     # kg/m³ or reduced
    end
    return (Tc, pc, ϱc)
end

end