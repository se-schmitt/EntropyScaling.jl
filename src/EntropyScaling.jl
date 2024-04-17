module EntropyScaling

# Load public modules
using LsqFit
using NLsolve
using Optim
using StatsBase
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
if !isdefined(Base,:get_extension)
    include("../ext/MicthermExt.jl")
    include("../ext/ClapeyronExt.jl")
end

# Export
export fit_entropy_scaling, call_entropy_scaling

# Export extension functions
# MicthermExt
export MicThermParam
function MicThermParam end
abstract type MicThermParamType end

end
