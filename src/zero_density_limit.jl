# Functions to calculate the scaled Chapman-Enskog (CE) transport properties

# Function to calculate scaled CE transport properties
function CE_scaled(m::NamedTuple, T::Vector{Float64}, prop::String; x=[], solute::Dict{Symbol,Float64}=Dict{Symbol,Float64}(), reduced::Bool=false, min_xdep::Bool=false)
    # Calculate LJ parameters
    σ_CE, ε_CE = calc_σε(m, prop; solute=solute, reduced=reduced)

    # Scaled transport property
    f =     prop == "vis" ? 5/16 :
            prop == "tcn" ? 75/64 :
            prop in ["dif","selfdif","mutdif"] ? 3/8 : error("'prop' must be 'vis', 'tcn', 'dif', 'selfdif', or 'mutdif'!")
            
    Ω(T) =  prop == "vis" ? Ω_22(T) :
            prop == "tcn" ? Ω_22(T) :
            prop in ["dif","selfdif","mutdif"] ? Ω_11(T) : error("'prop' must be 'vis', 'tcn', 'dif', 'selfdif', or 'mutdif'!")
    TdBdT_B(T,x) = isempty(x) ? ((T.*m.dBdTfun(T).+m.Bfun(T))./NA).^(2/3) :
                                ((T.*m.dBdTmixfun(T,x).+m.Bmixfun(T,x))./NA).^(2/3)

    index = prop == "vis" ? 1 :
            prop == "tcn" ? 2 :
            prop in ["dif","selfdif","mutdif"] ? 3 : error("'prop' must be 'vis', 'tcn', 'dif', 'selfdif', or 'mutdif'!")
    Y_CE⁺(T,x) = m.μ == 0.0 ?   f./((√(π)*σ_CE^2).*Ω(T./ε_CE.*kB)).*TdBdT_B(T,x) :
                                py"calc_Stockmayer_transport"(T,round(m.μ^2,digits=2))[index]

    # Calculate minimum of Y_CE⁺
    min_Y_CE⁺ = NaN
    
    TB = nlsolve(x -> m.Bfun(x[1]),[0.6*m.Tc]).zero[1]        # Boyle temperature
    try 
        opt = optimize(y -> Y_CE⁺(y[1],x),[TB],NewtonTrustRegion())
        min_Y_CE⁺ = opt.minimum
    catch e
        if isa(e,DomainError)
            @warn("DomainError in Y₀⁺! Used value at T = 0.6*T_Boyle as min(Y₀⁺).")
            min_Y_CE⁺ = Y_CE⁺(0.6*TB,x)
        else
            throw(e)
        end
    end

    return Y_CE⁺(T,x), min_Y_CE⁺
end

function calc_σε(m, prop; solute=Dict{Symbol,Float64}(), reduced=false)
    # Critical temperature and pressure of the LJ fluid from Stephan et al. (2019) [DOI: 10.1021/acs.jcim.9b00620]
    Tc_LJ = 1.321
    pc_LJ = 0.129

    if reduced
        # Take LJ parameters from model
        σ_CE = m.σ
        ε_CE = m.ε
        if !isempty(solute)
            ε_CE = sqrt(ε_CE * solute[:ε]) * solute[:ξ]
            σ_CE = (σ_CE + solute[:σ]) / 2
        end
    else
        # Apply correspondence principle to calculate LJ parameters
        ε_CE = kB.*m.Tc./Tc_LJ
        σ_CE = (pc_LJ./m.pc.*ε_CE).^(1/3)
        if prop in ["dif","mutdif","selfdif"] && !isempty(solute)
            ε_CE_sol = kB.*solute[:Tc]./Tc_LJ
            σ_CE_sol = (pc_LJ./solute[:pc].*ε_CE_sol).^(1/3)
            ε_CE = sqrt(ε_CE * ε_CE_sol)
            σ_CE = (σ_CE + σ_CE_sol) / 2
        end
    end

    # Calculate LJ parameters of mixture
    if all(length.([σ_CE, ε_CE]) .== 2)
        σ_CE = mean(σ_CE)
        ε_CE = geomean(ε_CE)
    elseif any(length.([σ_CE, ε_CE]) .≥ 2)
        error("Both σ_CE and ε_CE must be either scalar or 2-element vectors.")
    end

    return σ_CE, ε_CE
end

# Collision integrals from Kim and Monroe (2014) [DOI: 10.1016/j.jcp.2014.05.018]
# Ω₂₂
function Ω_22(T_red)
    A = -0.92032979
    BCi = [ [   2.3508044,      1.6330213       ],
            [   0.50110649,     -6.9795156e-1   ],
            [   -4.7193769e-1,  1.6096572e-1,   ],
            [   1.5806367e-1,   -2.2109440e-2   ],
            [   -2.6367184e-2,  1.7031434e-3    ],
            [   1.8120118e-3,   -0.56699986e-4    ]]
    out = A .+ sum([BCi[k][1] ./ T_red.^k .+ BCi[k][2] .* log.(T_red).^k for k in 1:6])
end

# Ω₁₁
function Ω_11(T_red)
    A = -1.1036729
    BCi = [ [   2.6431984,      1.6690746       ],
            [   0.0060432255,   -6.9145890e-1   ],
            [   -1.5158773e-1,  1.5502132e-1,   ],
            [   0.54237938e-1,  -2.0642189e-2   ],
            [   -0.90468682e-2, 1.5402077e-3    ],
            [   0.61742007e-3,  -0.49729535e-4  ]]
    out = A .+ sum([BCi[k][1] ./ T_red.^k .+ BCi[k][2] .* log.(T_red).^k for k in 1:6])
end