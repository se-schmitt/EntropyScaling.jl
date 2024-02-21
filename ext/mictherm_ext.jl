# Extension provides wrapper for functions call_entropy_scaling and fit_entropy_scaling to be coupled with the MicTherm package

using MATLAB

# MicTherm parameter struct
struct MicThermParam
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
    source::Union{Vector{String},Missing}
    path::String
end

"""
`MicThermParam(ParamDict::NamedTuple)`

Constructor for MicThermParam struct

- Mandatory input parameters: `:name`, `:eos`, `:path`, `:unit`
- Optional input parameters: `:polar`, `:assoc`, `:M`, `:m`, `:σ`, `:ε`, `:ξ`, `:Tcm`, `:pcm`, `:λr`, `:λa`, `:TStar`, `:vStar`, `:μ`, `:N_μ`, `:Q`, `:N_Q`, `:type_AB`, `:κ_AB`, `:ε_AB`, `:source`

Examplet for constructing a MicThermParam struct:
```julia
d = (name=["hexane", "heptane"], eos="PC_SAFT", path="/path/to/MicTherm", unit="SI", M=[86.18, 100.21], m=[3.0576,3.4831], σ=[3.7983,3.8049], ε=[236.77,238.4], source=["10.1021/ie0003887"])
model = MicThermParam(d)
```
"""
function MicThermParam(ParamDict::NamedTuple)
    # Check Dict
    mandatory_input = [:name, :eos, :path]
    for f in mandatory_input
        if !(f in keys(ParamDict))  error("Mandatory input parameter $f is missing.") end
    end

    fields = fieldnames(MicThermParam)
    return MicThermParam([f in keys(ParamDict) ? ParamDict[f] : missing for f in fields]...)
end

# Wrapper for the function `fit_entropy_scaling` to be used with MicTherm
function fit_entropy_scaling(   T::Vector{Float64}, 
                                ϱ::Vector{Float64}, 
                                Y::Vector{Float64}, 
                                prop::String; 
                                model::MicThermParam,
                                i_fit=[0,1,1,1,1])
    # Set MicTherm path
    eval_string("addpath(genpath('$(model.path)'))")
    mat"warning off"

    # Create function handles
    (sfun, Bfun, dBdTfun) = get_MicTherm_fun(model)

    # Calculate critical temperature and pressure
    (Tc, pc, ϱc) = calc_crit_MicTherm(model)

    m = ismissing(model.m) ? ones(length(model.name)) : model.m
    M = model.unit == "reduced" ? model.M : model.M./1e3

    return fit_entropy_scaling(T, ϱ, Y, prop; sfun=sfun, Bfun=Bfun[1], dBdTfun=dBdTfun[1], Tc=Tc[1], pc=pc[1], M=M, i_fit=i_fit, m_EOS=m)
end

# Wrapper for the function `call_entropy_scaling` to be used with MicTherm
function call_entropy_scaling(  T::Vector{Float64}, 
                                ϱ::Vector{Float64},
                                α_par::Vector{Vector{Float64}},
                                prop::String; 
                                x::Matrix{Float64}=ones(length(T),1),
                                model::MicThermParam)
    
    # Set MicTherm path
    eval_string("addpath(genpath('$(model.path)'))")
    mat"warning off"
    
    # Create function handles
    (sfun, Bfun, dBdTfun) = get_MicTherm_fun(model)

    # Calculate critical temperature and pressure
    (Tc, pc, ϱc) = calc_crit_MicTherm(model)

    m = ismissing(model.m) ? ones(length(model.name)) : model.m
    M = model.unit == "reduced" ? model.M : model.M./1e3

    return call_entropy_scaling(T, ϱ, α_par, prop; x=x, sfun=sfun, Bfun=Bfun, dBdTfun=dBdTfun, Tc=Tc, pc=pc, M=M, m_EOS=m)
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
            out *= string(set_param.(Ref(p),[:m,:σ,:ε],["chainlength_$j","sigma_$j","epsilon_$j"],i)...)
        elseif lowercase(p.eos) in ["scpa"]
            out *= string(set_param.(Ref(p),[:Tcm,:pcm],["Tcm_$j","pcm_$j"],i)...)
        elseif lowercase(p.eos) in ["pactplusb"]
            out *= string(set_param.(Ref(p),[:vStar,:TStar,:vStar,:TStar],["vStar_$j","TStar_$j","sigma_$j","epsilon_$j"],i)...)
        else
            error("Either (σ, ϵ, m) or (Tcm, pcm) or (vTStar, vStar) should be defined!")
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
    if N_comp == 2 & !ismissing(p.ξ)
        out *= ", 'Xsi_12 = $(ξ[1])'"
    end

    return out
end

# Function to provide function objectives from MicTherm
function get_MicTherm_fun(model)
    # Set conversion factors
    if model.unit == "reduced"
        Bconv = 1
    elseif model.unit == "SI"
        Bconv = 1e-3
    end

    # Configurational entropy
    sfun(T,ϱ,x) = 
    (   eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = sRes');");
        mat"[$names , ~ , $values] = IO_API( 'initialized', $ϱ, $T, [], $x );";                    
        return values[:,findfirst(names .== "sRes")[2]] )

    # Second virial coefficient
    Bfun = [T ->
        (   eval_string("IO_API( $(get_initialization_string(model; comp=i)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = B');");
            mat"[$names , ~ , $values] = IO_API( 'initialized', $(ones(size(T))), $T, [], $(ones(size(T))) );";                    
            return values[:,findfirst(names .== "B")[2] .* Bconv] 
        )
        for i in eachindex(model.name) ]
    
    # Derivative of second virial coefficient
    dBdTfun = [T -> 
        (   eval_string("IO_API( $(get_initialization_string(model)), 'calculationmode = API', 'APIMode = UserProperties', 'properties = B');");
            mat"[$names , ~ , $values_h] = IO_API( 'initialized', $(ones(size(T))), $(T.+h/2), [], $(ones(size(T))) );";
            mat"[$names , ~ , $values_l] = IO_API( 'initialized', $(ones(size(T))), $(T.-h/2), [], $(ones(size(T))) );";
            B_h = values_h[:,findfirst(name .== "B")[2]];
            B_l = values_l[:,findfirst(name .== "B")[2]];
            return (B_h .- B_l) ./ h .* Bconv
        )
        for i in eachindex(model.name) ]
    
    return sfun, Bfun, dBdTfun
end 

# Function to calculate critical point by MicTherm
function calc_crit_MicTherm(model)
    # Set conversion factors
    if model.unit == "reduced"
        pconv = 1
        ϱconv = 1
    elseif model.unit == "SI"
        pconv = 1e6
        ϱconv = 1.0./model.M
    end

    Tc = Float64[]
    pc = Float64[]
    ϱc = Float64[]
    for i in 1:length(model.name)
        eval_string("IO_API( $(get_initialization_string(model; comp=i)), 'calculationmode = VLE', 'dT = 100' )")
        mat"[$name_VLE , ~ , $values_VLE] = IO_API( 'initialized', [], [], [], [] )"
        mat"close all"
        push!(Tc,values_VLE[1,findfirst(name_VLE .== "T_VLE")[2]])          # K or reduced
        push!(pc,values_VLE[1,findfirst(name_VLE .== "p_VLE")[2]]*pconv)    # Pa or reduced
        push!(ϱc,values_VLE[1,findfirst(name_VLE .== "rho_L")[2]]*ϱconv)    # kg/m³ or reduced
    end
    return (Tc, pc, ϱc)
end