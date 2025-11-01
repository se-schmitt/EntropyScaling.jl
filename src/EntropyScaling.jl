module EntropyScaling

# Load public modules
using SimpleNonlinearSolve, Optimization, OptimizationOptimJL, StatsBase, ForwardDiff, DelimitedFiles
import LogExpFunctions #loaded by StatsBase.jl
import FillArrays #loaded by Optimization.jl
const Z1 = FillArrays.Fill(1,1)
# Definition of Constants
get_kBNAR() = (kB=1.380649e-23; NA=6.02214076e23; return (kB,NA,kB*NA))
const (kB, NA, R) = get_kBNAR()

const DB_PATH = normpath(Base.pkgdir(EntropyScaling),"database")

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
include("general/properties.jl")

# Utils
include("utils/data.jl")
include("utils/thermo.jl")
include("utils/misc.jl")
include("utils/database.jl")
include("utils/plots.jl")

# Models 
include("models/base.jl")
include("models/chapman_enskog.jl")
include("models/dilute_gas_correlations.jl")
include("models/framework.jl")
include("models/refprop_res.jl")

end
