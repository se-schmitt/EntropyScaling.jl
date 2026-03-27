using Test
using EntropyScaling
using Clapeyron
using CoolProp
using DelimitedFiles
using Unitful
using Symbolics

#TODO rename to test_[...].jl
include("framework.jl")
include("refprop_res.jl")
include("gces.jl")
include("chapman_enskog.jl")
include("unitful.jl")
include("symbolics.jl")
include("plotting.jl")