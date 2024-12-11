module EntropyScaling

# Load public modules
using SimpleNonlinearSolve, Optimization, StatsBase, ForwardDiff, DelimitedFiles

# Definition of Constants
get_kBNAR() = (kB=1.380649e-23; NA=6.02214076e23; return (kB,NA,kB*NA))
(kB, NA, R) = get_kBNAR()

const DB_PATH = normpath(Base.pkgdir(EntropyScaling),"data")

#equivalent to a' * b, but with general iterators
function _dot(a,b)
    res = zero(Base.promote_eltype(a,b))
    for i in 1:length(a)
        res += a[i]*b[i]
    end
    return res
end
# General
include("general/types.jl")
include("general/scalings.jl")
include("general/chapman_enskog.jl")

# Utils
include("utils/data.jl")
include("utils/thermo.jl")
include("utils/misc.jl")

# Models 
include("models/framework.jl")

# Extensions
if !isdefined(Base,:get_extension)
    include("../ext/MicthermExt.jl")
    include("../ext/ClapeyronExt.jl")
end

end
