# Main function to call entropy scaling
function call_entropy_scaling(  T::Vector{Float64}, 
                                ϱ::Vector{Float64},
                                α_par::Vector{Vector{Float64}},
                                prop::String; 
                                x::Matrix{Float64}=ones(length(T),1),
                                sfun::Function, Bfun::Vector, dBdTfun::Vector, 
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
    Y_CE⁺ = mix_Wilke(Y_CE⁺_all, x, M)
    min_Y_CE⁺ = mix_Wilke(repeat.(min_Y_CE⁺_all,length(T)), x, M)

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