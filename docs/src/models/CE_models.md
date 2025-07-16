# Chapman-Enskog Model

Chapman-Enskog model for transport properties at the zero-density limit based on the kinetic gas theory.

```math
\begin{aligned}
\eta_{\varrho \rightarrow 0}              &= \frac{5}{16} \sqrt{\frac{M k_{\rm B} T}{\pi N_{\rm A}}} \frac{1}{\sigma^2 \Omega^{(2,2)}} \\
\lambda_{\varrho \rightarrow 0}	          &= \frac{75}{64} k_{\rm B} \sqrt{\frac{R T}{M \pi}} \frac{1}{\sigma^2 \Omega^{(2,2)}}\\
D_{\varrho \rightarrow 0} \varrho^{\rm m} &= \frac{3}{8} \sqrt{\frac{M k_{\rm B} T}{\pi N_{\rm A}}} \frac{1}{\sigma^2 \Omega^{(1,1)}}
\end{aligned}
```

```@docs
EntropyScaling.ChapmanEnskogModel
EntropyScaling.Ω
```
