# Implementation of the Peng-Robinson (PR) equation of state for methane

# Read experiemntal data
header

# Parameters for methane
Tc = 190.4                                                      # K
pc = 4.6e6                                                      # Pa
ω = 0.011                                                       # -
M = 0.016043                                                    # kg mol⁻¹

# Constants
R = 8.31446261815324                                            # J K⁻¹ mol⁻¹

# General PR parameters
a = 0.45724*R^2*Tc^2/pc                                         # J m³ mol⁻²
b = 0.07780*R*Tc/pc                                             # m³ mol⁻¹
m = ω≤0.491 ? 0.37464+1.54226*ω+0.26992*ω^2 :                   # -
                0.379642+1.48503*ω+0.164423*ω^2+0.016666*ω^3
Δ₁ = 1+m*√(2)                                                   # -
Δ₂ = 1-m*√(2)                                                   # -
α(T) = (1+m*(1-√(T/Tc)))^2                                      # -
dαdT(T) = -m*(m*(1-√(T/Tc)+1))/(Tc*√(T/Tc))                     # -

# Pressure
pfun(T,ϱ) = (R*T)/(1/(ϱ/M)-b) - a*α(T)/(1/(ϱ/M)^2+2b/(ϱ/M)-b^2)             # Pa  

# Entropy
sfun(T,ϱ) = R*log(1-b*(ϱ/M)) -                                  # J K⁻¹ mol⁻¹
            a*dαdT(T)*log((Δ₁*b*(ϱ/M)+1)/(Δ₂*b*(ϱ/M)+1))/(b*m*(Δ₁-Δ₂))

# 2nd virial coefficient
Bfun(T) = b - a*α(T)/(R*T)                                      # m³ mol⁻¹

# Temperature deriviative of 2nd virial coefficient
dBdTfun(T) = a*α(T)/(R*T^2) - a*dαdT(T)/(R*T)                   # m³ mol⁻¹ K⁻¹

# Implicit calculation
using Roots

find_zeros.(pfun)



