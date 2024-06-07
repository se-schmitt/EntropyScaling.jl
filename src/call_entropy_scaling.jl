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
- `prop::String`: Transport property (`vis`, `tcn`, `dif`, `selfdif`, or `mutdif`)
- Keyword arguments:
    - `x::Matrix{Float64}`: Mole fractions of components (default: `x=ones(length(T),1)` -> only valid for pure substances)
    - `reduced::Bool`: Using LJ reduced units (default: `reduced=false`)
    - `difcomp::Int64`: Component for self-diffusion coefficient calculation (`difcomp ∈ {1,2}`) (default: `difcomp=0`)

Output:
- `Y::Vector{Float64}`: Transport property (η, λ, or D) in SI units (Pa s, W m⁻¹ K⁻¹, or m² s⁻¹)
"""
function call_entropy_scaling(  model::Dict{Symbol,Any},
                                T::Vector{Float64}, 
                                ϱ::Vector{Float64},
                                prop::String; 
                                x::Matrix{Float64}=ones(length(T),1),
                                reduced=false,
                                difcomp::Int64=0)
    
    # Check input
    check_input_call(T, ϱ, x, prop)
    m = get_model(model, prop, size(x,2))

    if reduced
        global (kB,NA,R) = (1.0,1.0,1.0)
    elseif (kB,NA,R) == (1.0,1.0,1.0)
        global (kB,NA,R) = get_kBNAR()
    end

    # Calculation of number density
    ϱN = ϱ ./ (x * m.M) .* NA                                 # [ϱN] = 1/m³

    # Overwrrite reference masses for diffusion coefficient
    M_ref = deepcopy(m.M)
    if prop in ["selfdif","mutdif"]
        M_dif = 2.0 ./ (1.0 ./ m.M[1] .+ 1.0 ./ m.M[2])
        if prop == "mutdif" 
            M_ref = repeat([M_dif],2)
        else prop == "selfdif"
            if difcomp ∉ [1,2]
                error("Specify component for self-diffusion coefficient calculation (`difcomp ∈ {1,2}`).")
            end
            M_ref[3-difcomp] = M_dif
        end
    end

    # Calculation of entropy
    s_conf = m.sfun(T,ϱ,x)
    m_mix = x * m.m_EOS
    s = -s_conf ./ R ./ m_mix

    # Calculation of scaled Chapman-Enskog (CE) transport properties and its minimum
    Y_CE⁺_all = []
    min_Y_CE⁺_all = []
    for i in 1:size(x,2)
        solute = Dict{Symbol,Float64}()
        if prop == "mutdif" || (prop == "selfdif" && i != difcomp)
            if !reduced
                solute[:Tc] = m.Tc[3-i]
                solute[:pc] = m.pc[3-i]
            else
                solute[:ε] = m.ε[3-i]
                solute[:σ] = m.σ[3-i]
                solute[:ξ] = m.ξ
            end
        end
        (Y_CE⁺_i, min_Y_CE⁺_i) = CE_scaled(split_m(m)[i], T, prop; solute=solute, reduced=reduced)
        push!(Y_CE⁺_all,Y_CE⁺_i)
        push!(min_Y_CE⁺_all,[min_Y_CE⁺_i])
    end
    if prop in ["vis","tcn"]
        Y_CE⁺ = mix_Wilke(Y_CE⁺_all, x, m.M)
        min_Y_CE⁺ = mix_Wilke(repeat.(min_Y_CE⁺_all,length(T)), x, m.M)
    elseif prop in ["dif","selfdif"]
        Y_CE⁺ = mix_Darken(Y_CE⁺_all, x)
        min_Y_CE⁺ = mix_Darken(repeat.(min_Y_CE⁺_all,length(T)), x)
    elseif prop in ["mutdif"]
        Y_CE⁺, _ = CE_scaled(m, T, prop; x=x, reduced=reduced)
        min_Y_CE⁺ = mix_Darken(repeat.(min_Y_CE⁺_all,length(T)), x)
    end

    # Calculation of transport property
    α_mix = mix_es_parameters(m.α, x)
    Y_fit = fun_es(s, α_mix; prop=prop)
    if prop in ["vis","dif","selfdif","mutdif"]
        Yˢ = exp.(Y_fit)
    else
        Yˢ = Y_fit
    end

    # Unscale transport property
    Y⁺ = Yˢ ./ (W(s)./Y_CE⁺ .+ (1.0 .- W(s))./min_Y_CE⁺)
    if prop == "vis"
        Y = Y⁺ .* ϱN.^(2/3).*sqrt.((x*m.M).*T*kB/NA) ./ (-s_conf/R).^(2/3)     
    elseif prop == "tcn"
        Y = Y⁺ .* ϱN.^(2/3).*kB.*sqrt.(R*T./(x*m.M)) ./ (-s_conf/R).^(2/3)
    elseif prop in ["dif","selfdif","mutdif"]
        Y = Y⁺ .* ϱN.^(-1/3).*sqrt.(R*T./(x*M_ref)) ./ (-s_conf/R).^(2/3)
    end

    return Y
end

# Function to check input
function check_input_call(T, ϱ, x, prop)
    # Check property value
    if !(prop in ["vis","tcn","dif","selfdif","mutdif"])
        error("Property must be 'vis', 'tcn', 'dif', 'selfdif', or 'mutdif'.")
    end

    if prop in ["selfdif","mutdif"] && !(size(x,2) == 2)
        error("Number of components must be 2 for properties 'selfdif' and 'mutdif'.")
    end

    # Check equal number of state points
    if !(length(T) == length(ϱ) == size(x,1))
        error("Length of T, ϱ, and x (# state points) must be equal.")
    end
end

# Function to check input model
function get_model(model_ori, prop, Ncomp; is_fit=false, reduced=false)
    model = deepcopy(model_ori)

    # Check mandatory keys
    required = [:sfun,:Bfun,:Tc,:pc,:M]
    if !all(in.(required,Ref(keys(model))))
        error("Fields missing in NamedTuple `model`: $(join(required[(!).(in.(required,Ref(keys(model))))], ", ")).")
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
    # dBdTmixfun
    if !haskey(model,:dBdTmixfun) && !is_fit
        Bfun = deepcopy(model[:Bmixfun])
        model[:dBdTmixfun] = (x,z) -> [ForwardDiff.derivative(y -> Bfun(y,z[i,:]'),x[i]) for i in eachindex(x)][length(x) == 1 ? 1 : (:)]
    end
    # EOS parameters
    if !haskey(model,:m_EOS)
        model[:m_EOS] = ones(length(model[:M]))
    end
    if reduced && any(!haskey.(Ref(model),[:ε,:σ]))
        error("Fields `ε` and `σ` must be specified in `model` for reduced units.")
    end
    if reduced && !haskey(model,:ξ)
        model[:ξ] = 1.0
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

# Mixing rule from Wilke (1950) [DOI: 10.1063/1.1747673] and Mason and Saxena (1958) [DOI: 10.1063/1.1724352]
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

# Mixing rule from Miller and Carman (1961) [DOI: 10.1039/TF9615702143]
function mix_Darken(D, x)
    Dinv = zeros(length(D[1]))
    for i in eachindex(D)
        Dinv += x[:,i]./D[i]
    end
    return 1.0 ./ Dinv
end

# Function to split the model
function split_m(m)
    Ncomp = length(m[:Tc])
    required = [:Bfun, :dBdTfun,:Tc,:pc]
    optional = [:Bmixfun, :dBdTmixfun, :ε, :σ]
    ms = Vector{NamedTuple}(undef,Ncomp)
    for i in 1:Ncomp
        id = Dict{Symbol,Any}()
        for k in required
            id[k] = m[k][i]
        end
        for k in optional
            if haskey(m,k)
                if k in [:ε,:σ]
                    id[k] = m[k][i]
                else
                    id[k] = m[k]
                end
            end
        end
        ms[i] = (;id...)
    end
    return ms
end
