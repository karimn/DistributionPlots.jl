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
include("dotlayout.jl")
include("geometry.jl")
include("recipes/slabinterval.jl")
# includes are added as later tasks create the files:
# include("recipes/dots.jl")
# include("recipes/lineribbon.jl")

export asdist
export point_interval
export slab_curve
export dot_layout
export slabinterval, slabinterval!, SlabInterval
export halfeye, halfeye!, eye, eye!, ccdfinterval, ccdfinterval!, cdfinterval, cdfinterval!,
       gradientinterval, gradientinterval!, histinterval, histinterval!,
       slab, slab!, interval, interval!, pointinterval, pointinterval!, spike, spike!,
       HalfEye, Eye, CcdfInterval, CdfInterval, GradientInterval, HistInterval,
       Slab, Interval, PointInterval, Spike

end # module
