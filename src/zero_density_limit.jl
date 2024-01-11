# Functions to calculate the scaled Chapman-Enskog (CE) transport properties

# Function to calculate scaled CE transport properties
function CE_scaled(T::Vector{Float64}, Tc::Float64, pc::Float64, prop::String, B::Function, dBdT::Function)
    # Critical temperature and pressure of the LJ fluid
    Tc_LJ = 1.321
    pc_LJ = 0.129

    # Apply correspondenve principle to calculate LJ parameters
    ε_CE = kB.*Tc./Tc_LJ
    σ_CE = (pc_LJ./pc.*ε_CE).^(1/3)

    # Scaled transport property
    f =     prop == "vis" ? 5/16 :
            prop == "tcn" ? 75/64 :
            prop == "dif" ? 3/8 : error("'prop' must be 'vis', 'tcn', or 'dif'!")
    Ω(T) =  prop == "vis" ? Ω_22(T) :
            prop == "tcn" ? Ω_22(T) :
            prop == "dif" ? Ω_11(T) : error("'prop' must be 'vis', 'tcn', or 'dif'!")
    Y_CE⁺(T) = f  / (√(π) * σ_CE^2 * Ω(T/ε_CE*kB)) * ((T[1]*dBdT(T)+B(T))/NA)^(2/3)

    # Calculate minimum of Y_CE⁺
    TB = nlsolve(x -> B(x[1]),[0.6*Tc]).zero[1]        # Boyle temperature
    min_Y_CE⁺ = NaN
    try 
        min_Y_CE⁺ = optimize(x -> Y_CE⁺(x[1]),[TB],NewtonTrustRegion()).minimum
    catch e
        if isa(e,DomainError)
            @warn("DomainError in Y₀⁺! Used value at T = 0.6*T_Boyle as min(Y₀⁺).")
            min_Y_CE⁺ = fun_η₀⁺(0.6*T_B)
        else
            throw(e)
        end
    end

    return Y_CE⁺.(T), min_Y_CE⁺
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