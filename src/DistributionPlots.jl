module DistributionPlots

using Statistics
using StatsBase
using Distributions
using KernelDensity
using Tables
using Makie

include("interface.jl")
include("density.jl")
include("intervals.jl")
include("slabs.jl")
# includes are added as later tasks create the files:
# include("dotlayout.jl")
# include("geometry.jl")
# include("recipes/slabinterval.jl")
# include("recipes/dots.jl")
# include("recipes/lineribbon.jl")

export asdist
export point_interval
export slab_curve

end # module
