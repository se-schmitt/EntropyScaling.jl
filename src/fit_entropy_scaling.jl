# Main function to fit entropy scaling parameters
function fit_entropy_scaling(   T::Vector{Float64}, 
                                ϱ::Vector{Float64}, 
                                Y::Vector{Float64}, 
                                prop::String; 
                                sfun::Function, Bfun::Function, dBdTfun::Function, 
                                Tc::Float64, pc::Float64, M::Float64,
                                i_fit=[0,1,1,1,1], m_EOS=1.0)

    # Calculation of entropy
    s_conf = sfun(T,ϱ,ones(length(T),1))
    s = -s_conf ./ R ./ m_EOS
    
    # Modified Rosenfeld scaling
    ϱN = ϱ ./ M .* NA                                               # [ϱN] = 1/m³
    if prop == "vis"
        Yʳ = @. Y / (ϱN^(2/3) * sqrt(M / NA * kB * T))
        Y⁺ = @. Yʳ * (-s_conf / R)^(2/3)  
    end
    if prop == "dif"
        Yʳ = @. Y * sqrt(M / (NA * kB * T)) * ϱN^(1/3)
        Y⁺ = @. Yʳ * (-s_conf / R)^(2/3)  
    end
    if prop == "tcn"
        Yʳ = @. Y / (ϱN^(2/3) * kB) * sqrt(M / (T * kB * NA)) 
        Y⁺ = @. Yʳ * (-s_conf / R)^(2/3)    
    end

    # Calculation of scaled Chapman-Enskog (CE) transport properties and its minimum
    (Y_CE⁺, min_Y_CE⁺) = CE_scaled(T, Tc, pc, prop, Bfun, dBdTfun)

    # CE-scaled transport properties
    Yˢ = (W(s)./Y_CE⁺ .+ (1.0 .- W(s))./min_Y_CE⁺) .* Y⁺

    # Fit component-specific parameters
    if prop in ["vis","dif"]    
        p0 = 0.
        Y_fit = log.(Yˢ)
    else                        
        p0 = 1.
        Y_fit = Yˢ
    end
    fit = curve_fit((x,par) -> fun_es(x,make_α(par,i_fit,p0);prop=prop), s, Y_fit, zeros(sum(i_fit)))
    α_par = make_α(fit.param,i_fit,p0)
    
    return α_par, Yˢ, s
end

# Function to create fit vector 
make_α(p,i,p0) = setindex!(vcat(p0,zeros(4)),p,i .== 1)