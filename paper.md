---
title: 'EntropyScaling.jl: Consistent thermodynmic modelling of transport properties in Julia'
tags:
  - Julia
  - thermodynamics
  - transport properties
  - viscosity
  - thermal conductivity
  - diffusion
  - equations of state (EOS)
authors:
  - name: Sebastian Schmitt
    orcid: 0000-0001-7909-2945
    affiliation: 1
affiliations:
 - name: Laboratory of Engineering Thermodynamics (LTD), RPTU Kaiserslautern
   index: 1
date: 01 April 2026
bibliography: paper.bib
---

# Summary

Tranport properties of fluids -- like the viscosity, the thermal conductivity, and diffusion coefficients -- are crucial for the design of variuos industrial processes.
For example, knowledge of the viscostity is required for accurate fluid dynamics simulations, the thermal conductivity is a central property for modeling heat transfer, and models for diffusion coefficients are used for the rate-based design of separation processes.
Accurate models of the transport properties enable the optimization of these processes in combination with physical simulations.
Conventional methods to model transport properties are often simple correlations that are limited to the scope of the data used in the fitting process.
In contrast, entropy scaling provides a very physical method for modeling transport properties that utilizes the fact that the transport properties, if appropriatly reduced, form a univariate of the residual entropy of the fluid.
This observation can be used to robustly model transport properties in all fluid states (gas, liquid, supercritical) by correlating this univariate function whereas the residual entropy is calculated from fundamental equations of state (EOS).
Fundamental EOS are formulations of the (residual) Helmholtz energy which enable the calculation of all (static) thermodynamic properties.
They are poweful models often used in thermodynamics.
In recent years, several models were developed to model transport properties based on entropy scaling.

The Julia package `EntropyScaling.jl` provides easily accessible implementations of entropy scaling models and thus enables the prediction of transport properties in different states in a physical manner.
It includes predefined models for some fluids, group-contribution models that enable the prediction of transport properties of components without experimental data, as well as methods to fit component-specific models with only a minimal amount of experimental data.
Additionally, `EntropyScaling.jl` includes an implementation of the Chapman-Enskog theory to calculate the transport properties of dilute gases.
The rich thermodynamics library `Clapeyron.jl` [@walker_clapeyronjl_2022] is used as a backend for the EOS calculations.

# Statement of need

`EntropyScaling.jl` is a thermodynamic modeling package for transport properties like the viscosity, thermal conductivity, or diffusion coefficients based on physically sound methods.
Among those, entropy scaling, which is also giving the package its name, enables modeling transport properties in a wide range of different states while only requiring a minimal set of experimental data by exploting the predictive capabilities of equations of state.
Besides entropy scaling, transport property models for gases based on the Chapman-Enskog theory are implemented.
Through `EntropyScaling.jl`, those methods for modeling transport properties are made available for applications in engineering tasks, e.g. separation processes or heat transfer.
Through Julia's excellent extensibility, the package can easily be coupled with the wider modeling ecosystem.
Moreover, using the strong interoperability with other programming languages (in particular Python), the methods from `EntropySclaing.jl` can also be used in a large number of different applications.

# Key features

## Scope

- list properties that can be modeled
- list "extras" like fitting functionality, plotting

## Available Models

- list all models with brief descripion

## Examples

- extract two example (incl. a nice plot, see @schmitt_entropy_2024 or @schmitt_entropy_2025 for example plots -> would be nice to have the code for a published code (take one plot out of every paper, viscosity from 2024, MS+Fick from 2025))

# State of the field                                                                                                                  

There are few packages that implement methods for modeling transport properties.
CoolProp [@bell_pure_2014] implements highly accurate correlations for the viscosity and thermal conductivity of few fluids for which a large number of experimental exist.
FeOS [@rehner_feos_2023] is Rust package with a Python frontend for thermodynamic calculations based on molecular equations of state like PC-SAFT EOS [@gross_perturbed_2001].
Additionally, it provides methods for calculating the viscosity, thermal conductivity, and self-diffusion coefficients based on entropy scaling.
Thermo [@bell_thermo_2016] is a general thermodynamic dynamic library written in Python that implements, besides other models for static thermodynamic property prediction, some simple correlations for the viscosity and thermal conductivity as well as the Joback group-contribution model [@joback_estimation_1987].
In Clapeyron.jl [@walker_clapeyronjl_2022], only the Joback group contribution is available for predicting the viscosity.

# Software design

- extension of clapeyron / tight integration + similar syntax
- easy to implement new models -> ES models only need few new functions 

# Research impact statement

- used in publications [@schmitt_entropy_2025]
- integrated in MLThermoProperties package (link)

# AI usage disclosure

The authors used generative AI tools in the preparation of both this manuscript and the EntropyScaling.jl package.

For the manuscript, AI tools were used to improve the initial draft. All AI-generated text was subsequently reviewed, edited, and integrated by the human authors.

For the software, AI tools provided supporting assistance only: code review, minor tasks, and testing. The design and the main body of the code were written by the authors, and any AI-generated contribution was inspected and tested before inclusion in the repository.

The authors are able to trace and justify all AI-assisted output and are responsible for the final content and its accuracy.

# Acknowledgements

- write sth general

# References