module EntropyScaling

# Load public modules
using LsqFit
using NLsolve
using Optim
using Statistics
using ForwardDiff

# Definition of Constants
get_kBNAR() = (kB=1.380649e-23; NA=6.02214076e23; return (kB,NA,kB*NA))
(kB, NA, R) = get_kBNAR()

# Include files
include("fit_entropy_scaling.jl")
include("call_entropy_scaling.jl")
include("correlation_fun.jl")
include("zero_density_limit.jl")

# Include extensions
if :MATLAB in names(Main, imported=true)
    include("../ext/mictherm_ext.jl")
end

# Export
export fit_entropy_scaling, call_entropy_scaling

end
