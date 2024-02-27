# Main function to call entropy scaling
"""
`call_entropy_scaling(model, T, ϱ, prop; x, reduced)`

Function to calculate transport properties using entropy scaling.

---

Input:
- `model::Dict{Symbol,Any}`: Dict with model parameters with keys:
    - `sfun::Function`: Function to calculate entropy in SI units (J K⁻¹ mol⁻¹) (`sfun(T,ϱ,x)`)
    - `Bfun::Vector`: Functions to calculate 2nd virial coefficient of the pure components in m³ mol⁻¹ (`Bfun(T)`)
    - `dBdTfun::Vector`: Functions to calculate temperature derivative of 2nd virial coefficient of the pure components in m³ mol⁻¹ K⁻¹ (`dBdTfun(T)`) (optional, calculated by automatic differentiation from `Bfun` otherwise)
    - `Tc::Vector{Float64}`: Critical temperatures of the pure components in K
    - `pc::Vector{Float64}`: Critical pressures of the pure components in Pa
    - `M::Vector{Float64}`: Molar masses of the pure components in kg mol⁻¹
    - `m_EOS::Vector{Float64}`: Segment numbers of the pure components (if not specified, `m_EOS=ones(length(α_par))`)
    - `α::Vector{Vector{Float64}}`: Component-specific parameters α₀ - α₄
- `T::Vector{Float64}`: Temperature in K
- `ϱ::Vector{Float64}`: Density in kg/m³
- `prop::String`: Transport property (`vis`, `tcn`, or `dif`)
- Keyword arguments:
    - `x::Matrix{Float64}`: Mole fractions of components (default: `x=ones(length(T),1)` -> only valid for pure substances)
    - `reduced::Bool`: Using LJ reduced units (default: `reduced=false`)

Output:
- `Y::Vector{Float64}`: Transport property (η, λ, or D) in SI units (Pa s, W m⁻¹ K⁻¹, or m² s⁻¹)
"""
function call_entropy_scaling(  model::Dict{Symbol,Any},
                                T::Vector{Float64}, 
                                ϱ::Vector{Float64},
                                prop::String; 
                                x::Matrix{Float64}=ones(length(T),1),
                                reduced=false)
    
    # Check input
    check_input_call(T, ϱ, x, prop)
    m = get_model(model, prop, size(x,2))

    if reduced
        global (kB,NA,R) = (1.0,1.0,1.0)
    elseif (kB,NA,R) == (1.0,1.0,1.0)
        global (kB,NA,R) = get_kBNAR()
    end

    # Calculation of entropy
    s_conf = m.sfun(T,ϱ,x)
    m_mix = x * m.m_EOS
    s = -s_conf ./ R ./ m_mix

    # Calculation of scaled Chapman-Enskog (CE) transport properties and its minimum
    Y_CE⁺_all = []
    min_Y_CE⁺_all = []
    for i in 1:size(x,2)
        (Y_CE⁺_i, min_Y_CE⁺_i) = CE_scaled(T, m.Tc[i], m.pc[i], prop, m.Bfun[i], m.dBdTfun[i])
        push!(Y_CE⁺_all,Y_CE⁺_i)
        push!(min_Y_CE⁺_all,[min_Y_CE⁺_i])
    end
    Y_CE⁺ = mix_Wilke(Y_CE⁺_all, x, m.M)
    min_Y_CE⁺ = mix_Wilke(repeat.(min_Y_CE⁺_all,length(T)), x, m.M)

    # Calculate of transport property
    α_mix = mix_es_parameters(m.α, x)
    Y_fit = fun_es(s, α_mix; prop=prop)
    if prop in ["vis","dif"]
        Yˢ = exp.(Y_fit)
    else
        Yˢ = Y_fit
    end

    # Unscale transport property
    Y⁺ = Yˢ ./ (W(s)./Y_CE⁺ .+ (1.0 .- W(s))./min_Y_CE⁺)
    ϱN = ϱ ./ (x * m.M) .* NA                                 # [ϱN] = 1/m³
    M_mix = x * m.M
    if prop == "vis"
        Y = Y⁺ .* ϱN.^(2/3).*sqrt.((x*m.M).*T*kB/NA) ./ (-s_conf/R).^(2/3)     
    elseif prop == "tcn"
        Y = Y⁺ .* ϱN.^(2/3).*kB.*sqrt.(R*T./(x*m.M)) ./ (-s_conf/R).^(2/3)
    elseif prop == "dif"
        Y = Y⁺ .* ϱN.^(-1/3).*sqrt.(R*T./(x*m.M)) ./ (-s_conf/R).^(2/3)
    end

    return Y
end

# Function to check input
function check_input_call(T, ϱ, x, prop)
    # Check property value
    if !(prop in ["vis","tcn","dif"])
        error("Property must be 'vis', 'tcn', or 'dif'.")
    end

    # Check equal number of state points
    if !(length(T) == length(ϱ) == size(x,1))
        error("Length of T, ϱ, and x (# state points) must be equal.")
    end
end

# Function to check input model
function get_model(model_ori, prop, Ncomp; is_fit=false)
    model = deepcopy(model_ori)

    # Check mandatory keys
    required = [:sfun,:Bfun,:Tc,:pc,:M]
    if !all(in.(required,Ref(keys(model))))
        error("Fields missing in NamedTuple `model`: $(join(required[(!).(in.(required,keys(model)))], ", ")).")
    end
    if !haskey(model,:α) && !is_fit
        sym = prop == "vis" ? :α_η : prop == "tcn" ? :α_λ : prop == "dif" ? :α_D : []
        if haskey(model,sym)
            model[:α] = model[sym]
        else
            error("Field `$(sym)` missing in NamedTuple `model`.")
        end
    end

    # Check and set optional keys
    # dBdTfun
    if !haskey(model,:dBdTfun)
        Bfun = deepcopy(model[:Bfun])
        if isa(Bfun,Function)
            model[:dBdTfun] = x -> ForwardDiff.derivative.(Bfun,x)
        else
            model[:dBdTfun] = [x -> ForwardDiff.derivative.(f,x) for f in Bfun]
        end
    end
    # m_EOS
    if !haskey(model,:m_EOS)
        model[:m_EOS] = ones(length(model[:M]))
    end

    # Check and set vector keys
    if !is_fit
        vector_keys = [:Bfun,:dBdTfun,:Tc,:pc,:M,:m_EOS,:α]
        if Ncomp == 1
            for k in vector_keys
                if any(isa.(Ref(model[k]),[Function,Number]))
                    model[k] = [model[k]]
                elseif k == :α && isa(model[k],Vector{Float64})
                    model[k] = [model[k]]
                end
            end
        end
        if (!).(all(length.([model[k] for k in vector_keys]) .== Ncomp))
            error("Lengths of fields in Dict `model` must be equal.")
        end
    end

    return (; model...)
end

# Function to calculate paramters of mixture
function mix_es_parameters(α_par::Vector{Vector{Float64}}, x::Matrix{Float64})
    # Mixing of parameters
    α_mix = Vector{Vector{Float64}}(undef,0)
    for i in 1:5
        push!(α_mix,x * getindex.(α_par,i))
    end
    
    return α_mix
end

# Mixing rule of Wilke 
function mix_Wilke(Y, x, M)
    Y_mix = zeros(length(Y[1]))
    for i in eachindex(Y)
        xΦ = zeros(length(Y[1]))
        for j in eachindex(Y)
            xΦ += x[:,j] .* (1.0 .+ (Y[i]./Y[j]).^(1/2) .* (M[j]./M[i]).^(1/4)).^2 ./ (8.0.*(1.0 .+ M[i]./M[j])).^(1/2)
        end
        Y_mix += x[:,i].*Y[i] ./ xΦ
    end
    return Y_mix
end