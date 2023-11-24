# Peng-Robinson equation of state
# The function PR() provides functions to calculate
#   - the pressure p_PR(T,ϱ)
#   - the entropy s_PR(T,ϱ)
#   - the 2nd virial coefficient B_PR(T)
#   - the temperature derivative of the 2nd virial coefficient dBdT_PR(T)
# of the Peng-Robinson (PR) equation of state for a given temperature T, density ϱ, critical temperature Tc, critical pressure pc, and acentric factor ω.
# Input variables are
#   - the critical temperature Tc with [Tc] = K
#   - the critical pressure pc with [pc] = Pa
#   - and the acentric factor ω with [ω] = 1
#   - the molar mass M with [M] = kg mol⁻¹
function PR(Tc, pc, ω, M)
    # Constants
    R = 8.31446261815324                                            # J K⁻¹ mol⁻¹

    # General PR parameters
    a = 0.45724*R^2*Tc^2/pc                                         # J m³ mol⁻²
    b = 0.07780*R*Tc/pc                                             # m³ mol⁻¹
    m = ω≤0.491 ?   0.37464+1.54226*ω+0.26992*ω^2 :                 # -
                    0.379642+1.48503*ω+0.164423*ω^2+0.016666*ω^3
    Δ₁ = 1+√(2)                                                     # -
    Δ₂ = 1-√(2)                                                     # -
    α(T) = (1+m*(1-√(T/Tc)))^2                                      # -
    dαdT(T) = -m*(m*(1-√(T/Tc)+1))/(Tc*√(T/Tc))                     # -

    # Pressure
    p_PR(T,ϱ,x) = @. (R*T)/(1/(ϱ/M)-b) -                         # Pa
                     a*α(T)/(1/(ϱ/M)^2+2b/(ϱ/M)-b^2)

    # Entropy
    s_PR(T,ϱ,x) = @. R*log(1-b*(ϱ/M)) +                             # J K⁻¹ mol⁻¹
                     a*dαdT(T)*log((Δ₁*b*(ϱ/M)+1.0)/(Δ₂*b*(ϱ/M)+1))/(b*m*(Δ₁-Δ₂))

    # 2nd virial coefficient
    B_PR(T) = @. b - a*α(T)/(R*T)                                   # m³ mol⁻¹

    # Temperature deriviative of 2nd virial coefficient
    dBdT_PR(T) = @. a*α(T)/(R*T^2) - a*dαdT(T)/(R*T)                # m³ mol⁻¹ K⁻¹

    return p_PR, s_PR, B_PR, dBdT_PR
end