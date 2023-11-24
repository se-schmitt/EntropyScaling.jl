# EntropyScaling

[![Build Status](https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/se-schmitt/EntropyScaling.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This is an implementation of a general entropy scaling framework intorduced in

S. Schmitt, H. Hasse, S. Stephan, Entropy Scaling Framework for Transport Properties using Molecular-based Equations of State, *Molecular Liquids* (2023) submitted.

## Installation

To register the module locally, type 
```julia
Pkg> add https://github.com/se-schmitt/EntropyScaling.jl
```
in package mode (type `]` to enter to Pkg mode in the REPL).

Then, the module can be loaded by
```julia
using EntropyScaling
```

## Usage

The module provides to main functions: `fit_entropy_scaling` and `call_entropy_scaling`.

### Function `fit_entropy_scaling`

**Input**:
- `T::Vector{Float64}`: temperature $T$ vector ($[T] = {\rm K}$) 
- `ϱ::Vector{Float64}`: density $\rho$ vector ($[\rho] = {\rm kg\,m^{-3}}$) 
- `Y::Vector{Float64}`: transport property $Y$ vector, can be
  - viscosity $\eta$ ($[\eta] = {\rm Pa\,s}$)
  - thermal conductivity $\lambda$ ($[\lambda] = {\rm W\,m^{-1}\,K^{-1}}$)
  - self-diffusion coefficient $D$ ($[D] = {\rm m^2\,s^{-1}}$)
- `prop::String`: property string, either `vis` for Viscosity, `tcn` for thermal conductivity, or `D` for self-diffusion coefficient
- `sfun::Function`: function to calculate the entropy $s$ ($[s] = {\rm J\,K^{-1}\,mol^{-1}$) of the form `sfun(T,ϱ,x)` (vectorized)
- `Bfun::Function`: function to calculate the 2nd virial coefficient $B$ ($[B] = {\rm m^3\,mol^{-1}$) of the form `Bfun(T)` (vectorized)
- `dBdTfun::Function`: function to calculate the temperature derivative of the 2nd virial coefficient ${\rm d} B/{\rm d} T$ ($[{\rm d} B/{\rm d} T] = {\rm m^3\,mol^{-1}\,K^{-1}$) of the form `dBdTfun(T)` (vectorized)
- `Tc::Float64`: critical temperature $ T_{\rm c} $ ($[T_{\rm c}] = {\rm K}$)
- `pc::Float64`: critical pressure $ p_{\rm c} $ ($[p_{\rm c}] = {\rm Pa}$)
- `M::Float64`: molar mass $M$ ($[M] = {\rm kg\,mol^{-1}}$)
- `i_fit`: vector determining which parameters of the entropy scaling model should be fitted $(\alpha_{0,i}, \alpha_{{\rm ln},i}, \alpha_{1,i}, \alpha_{2,i}, \alpha_{3,i})$ (if not specified: `i_fit=[0,1,1,1,1]`)
- `m_EOS`: segment number of the applied EOS (if not specified, `m_EOS=1.0`) 

**Output**:
- `α_par`: Fitted component specific parameters $(\alpha_{0,i}, \alpha_{{\rm ln},i}, \alpha_{1,i}, \alpha_{2,i}, \alpha_{3,i})$
- `Yˢ`: CE-scaled transport property (dimensionless), either $\widehat{\eta}^{+}$, $\widehat{\lambda}^{+}$, or $\widehat{D}^{+}$
- `s`: reduced configurational entropy `\tilde{s}_{\rm conf}`

### Function `call_entropy_scaling`



### Fit

### Units

Units applied in the module:

| Property | Symbol | Unit |
| --- | --- | --- |
| Temperature | $T$ | K |
| Pressure | $p$ | Pa |
| Density | $rho$ | kg/m³ |