# Main function to fit entropy scaling parameters
"""
`fit_entropy_scaling(T, ϱ, Y, prop; sfun, Bfun, dBdTfun, Tc, pc, M, i_fit, m_EOS, solute)`

Function to fit component-specific entropy scaling parameters for transport properties.

---

Input:
- `T::Vector{Float64}`: Temperature in K
- `ϱ::Vector{Float64}`: Density in kg/m³
- `Y::Vector{Float64}`: Transport property (η, λ, or D) in SI units (Pa s, W m⁻¹ K⁻¹, or m² s⁻¹)
- `prop::String`: Transport property (`vis`, `tcn`, or `dif`)
- Keyword arguments:
    - `sfun::Function`: Function to calculate entropy in SI units (J K⁻¹ mol⁻¹) (`sfun(T,ϱ,x)`)
    - `Bfun::Function`: Function to calculate 2nd virial coefficient in m³ mol⁻¹ (`Bfun(T)`)
    - `dBdTfun::Function`: Function to calculate temperature derivative of 2nd virial coefficient in m³ mol⁻¹ K⁻¹ (`dBdTfun(T)`) (optional, calculated by automatic differentiation from `Bfun` otherwise)
    - `Tc::Float64`: Critical temperature in K
    - `pc::Float64`: Critical pressure in Pa
    - `M::Float64`: Molar mass in kg mol⁻¹
    - `i_fit::Vector{Int64}`: Vector to specify which parameters should be fitted (default: `i_fit=[0,1,1,1,1]`) [length: 5]
    - `m_EOS::Float64`: segment number of the applied EOS (if not specified, `m_EOS=1.0`) 
    - `solute::Dict{Symbol,Float64}`: Dict with solute properties `:M` (molar mass), `:Tc` (critical temperature), and `:pc` (critical pressure) (optional, only relevant for `prop="dif"`)

Output:
- `α_par::Vector{Float64}`: Fitted component specific parameters α₀ - α₄
- `Yˢ::Vector{Float64}`: CE-scaled transport property (dimensionless), either ηˢ, λˢ, or Dˢ
- `s::Vector{Float64}`: Reduced configurational entropy (dimensionless)
"""
function fit_entropy_scaling(   T::Vector{Float64}, 
                                ϱ::Vector{Float64}, 
                                Y::Vector{Float64}, 
                                prop::String; 
                                sfun::Function, 
                                Bfun::Function, 
                                dBdTfun::Function=(x -> @. ForwardDiff.derivative(Bfun,x)), 
                                Tc::Float64, pc::Float64, M::Float64,
                                i_fit=[0,1,1,1,1], m_EOS=1.0,
                                solute::Dict{Symbol,Float64}=Dict{Symbol,Float64}())

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
        if !isempty(solute)
            M = 2/(1/M + 1/solute[:M])
        end
        Yʳ = @. Y * sqrt(M / (NA * kB * T)) * ϱN^(1/3)
        Y⁺ = @. Yʳ * (-s_conf / R)^(2/3)  
    end
    if prop == "tcn"
        Yʳ = @. Y / (ϱN^(2/3) * kB) * sqrt(M / (T * kB * NA)) 
        Y⁺ = @. Yʳ * (-s_conf / R)^(2/3)    
    end

    # Calculation of scaled Chapman-Enskog (CE) transport properties and its minimum
    (Y_CE⁺, min_Y_CE⁺) = CE_scaled(T, Tc, pc, prop, Bfun, dBdTfun; solute=solute)

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