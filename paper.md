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



# State of the field                                                                                                                  



# Software design



# Research impact statement



#-------------------------------------------------------------------------------------------

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

#-------------------------------------------------------------------------------------------

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References