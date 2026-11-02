module EntropyScaling

# Load public modules
using Optim, LsqFit, StatsBase, DelimitedFiles, LinearAlgebra
using ADTypes: AutoForwardDiff
using ForwardDiff: ForwardDiff
using FillArrays: FillArrays
using Clapeyron: Clapeyron
const CL = Clapeyron
const Z1 = FillArrays.Fill(1, 1)

# Definition of constants
get_kBNAR() = (kB=1.380649e-23; NA=6.02214076e23; return (kB, NA, kB*NA))
const (kB, NA, R) = get_kBNAR()

# k_B, N_A and the factors from the input units of σ (Å) and Mw (g mol⁻¹) to SI.
# Clapeyron potential models work in Lennard-Jones units, where all of them are 1.
const SI_UNITS = (kB=kB, NA=NA, σ=1e-10, Mw=1e-3)
const LJ_UNITS = (kB=1.0, NA=1.0, σ=1.0, Mw=1.0)
unit_constants(eos) = SI_UNITS
unit_constants(::CL.PotentialModel) = LJ_UNITS

const DB_PATH = normpath(Base.pkgdir(EntropyScaling), "database")

# equivalent to a' * b, but with general iterators
function _dot(a, b)
    res = zero(Base.promote_eltype(a, b))
    for i in 1:length(a)
        res += a[i] * b[i]
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
include("utils/split_models.jl")

# Models
include("models/base.jl")
include("models/utils.jl")
include("models/chapman_enskog.jl")
include("models/dilute_gas_correlations.jl")
include("models/framework.jl")
include("models/refpropres.jl")
include("models/refpropres2025.jl")
include("models/gces.jl")

end
