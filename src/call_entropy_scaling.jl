# Main function to call entropy scaling
"""
`call_entropy_scaling(T, ϱ, α_par, prop; x, sfun, Bfun, dBdTfun, Tc, pc, M, m_EOS)`

Function to calculate transport properties using entropy scaling.

---

Input:
- `T::Vector{Float64}`: Temperature in K
- `ϱ::Vector{Float64}`: Density in kg/m³
- `α_par::Vector{Vector{Float64}}`: Component-specific parameters α₀ - α₄
- `prop::String`: Transport property (`vis`, `tcn`, `selfdif` (`N ≤ 2`), or `mutdif` (`N ≤ 2`))
- Keyword arguments:
    - `x::Matrix{Float64}`: Mole fractions of components (default: `x=ones(length(T),1)` -> only valid for pure substances)
    - `sfun::Function`: Function to calculate entropy in SI units (J K⁻¹ mol⁻¹) (`sfun(T,ϱ,x)`)
    - `Bfun::Vector`: Functions to calculate 2nd virial coefficient of the pure components in m³ mol⁻¹ (`Bfun(T)`)
    - `dBdTfun::Vector`: Functions to calculate temperature derivative of 2nd virial coefficient of the pure components in m³ mol⁻¹ K⁻¹ (`dBdTfun(T)`) (optional, calculated by automatic differentiation from `Bfun` otherwise)
    - `Tc::Vector{Float64}`: Critical temperatures of the pure components in K
    - `pc::Vector{Float64}`: Critical pressures of the pure components in Pa
    - `M::Vector{Float64}`: Molar masses of the pure components in kg mol⁻¹
    - `m_EOS::Vector{Float64}`: Segment numbers of the pure components (if not specified, `m_EOS=ones(length(α_par))`)

Output:
- `Y::Vector{Float64}`: Transport property (η, λ, or D) in SI units (Pa s, W m⁻¹ K⁻¹, or m² s⁻¹)
"""
function call_entropy_scaling(  T::Vector{Float64}, 
                                ϱ::Vector{Float64},
                                α_par::Vector{Vector{Float64}},
                                prop::String; 
                                x::Matrix{Float64}=ones(length(T),1),
                                sfun::Function, 
                                Bfun::Vector, 
                                dBdTfun::Vector=[(x -> @. ForwardDiff.derivative(f,x)) for f in Bfun], 
                                Tc::Vector{Float64}, pc::Vector{Float64}, M::Vector{Float64},
                                m_EOS=ones(length(α_par)))
    
    # Check input
    check_input_call(T, ϱ, x, α_par, prop, Bfun, dBdTfun, Tc, pc, M, m_EOS)

    # Calculation of entropy
    s_conf = sfun(T,ϱ,x)
    m_mix = x * m_EOS
    s = -s_conf ./ R ./ m_mix

    # Calculation of scaled Chapman-Enskog (CE) transport properties and its minimum
    Y_CE⁺_all = []
    min_Y_CE⁺_all = []
    for i in size(x,2)
        (Y_CE⁺_i, min_Y_CE⁺_i) = CE_scaled(T, Tc[i], pc[i], prop, Bfun[i], dBdTfun[i])
        push!(Y_CE⁺_all,Y_CE⁺_i)
        push!(min_Y_CE⁺_all,[min_Y_CE⁺_i])
    end
    if prop in ["vis","tcn"]
        Y_CE⁺ = mix_Wilke(Y_CE⁺_all, x, M)
        min_Y_CE⁺ = mix_Wilke(repeat.(min_Y_CE⁺_all,length(T)), x, M)
    else
        Y_CE⁺ = mix_Darken(Y_CE⁺_all, x)
        min_Y_CE⁺ = mix_Darken(repeat.(min_Y_CE⁺_all,length(T)), x)
    end

    # Calculate of transport property
    α_mix = mix_es_parameters(α_par, x)
    Y_fit = fun_es(s, α_mix; prop=prop)
    if prop in ["vis","dif"]
        Yˢ = exp.(Y_fit)
    else
        Yˢ = Y_fit
    end

    # Unscale transport property
    Y⁺ = Yˢ ./ (W(s)./Y_CE⁺ .+ (1.0 .- W(s))./min_Y_CE⁺)
    ϱN = ϱ ./ (x * M) .* NA                                 # [ϱN] = 1/m³
    if prop == "vis"
        Y = @. Y⁺ * ϱN^(2/3)*sqrt((x*M)*T*kB/NA) / (-s_conf/R)^(2/3)     
    elseif prop == "tcn"
        Y = @. Y⁺ * ϱN^(2/3)*kB*sqrt(R*T/(x*M)) / (-s_conf/R)^(2/3)
    elseif prop == "dif"
        Y = @. Y⁺ * ϱN^(-1/3)*sqrt(R*T/(x*M))  / (-s_conf/R)^(2/3)
    end

    return Y
end

# Function to check input
function check_input_call(T, ϱ, x, α_par, prop, Bfun, dBdTfun, Tc, pc, M, m_EOS)
    # Check property value
    if !(prop in ["vis","tcn","dif"])
        error("Property must be 'vis', 'tcn', or 'dif'.")
    end

    # Check equal number of state points
    if !(length(T) == length(ϱ) == size(x,1))
        error("Length of T, ϱ, and x (# state points) must be equal.")
    end
    
    # Check equal number of components
    if !(size(x,2) == length(α_par) == length(Bfun) == length(dBdTfun) == length(Tc) == length(pc) == length(M) == length(m_EOS))
        error("Number of components in  must be equal.")
    end
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