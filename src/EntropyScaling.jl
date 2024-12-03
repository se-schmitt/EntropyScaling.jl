module EntropyScaling

# Load public modules
using LsqFit
using NonlinearSolve
using Optimization
using StatsBase
using ForwardDiff

# Definition of Constants
get_kBNAR() = (kB=1.380649e-23; NA=6.02214076e23; return (kB,NA,kB*NA))
(kB, NA, R) = get_kBNAR()

# General
include("general/types.jl")
include("general/scalings.jl")
include("general/chapman_enskog.jl")

# Utils
include("utils/data.jl")
include("utils/thermo.jl")

# Models 
include("models/framework.jl")

# Extensions
if !isdefined(Base,:get_extension)
    include("../ext/MicthermExt.jl")
    include("../ext/ClapeyronExt.jl")
end

end
