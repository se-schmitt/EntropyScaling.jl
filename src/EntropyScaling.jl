module EntropyScaling

# Load public modules
using LsqFit
using NLsolve
using Optim
using Statistics

# Definition of Constants
kB = 1.380649e-23
NA = 6.02214076e23
R = kB*NA

# Include files
include("fit_entropy_scaling.jl")
include("call_entropy_scaling.jl")
include("correlation_fun.jl")
include("zero_density_limit.jl")

# Export
export fit_entropy_scaling, call_entropy_scaling

end
