module DistributionPlots

using Statistics
using StatsBase
using Distributions
using KernelDensity
using Tables
using Makie

"""
    dimticks(x, dim = 1)

Makie tick specification `(positions, labels)` for a labelled input. Implemented by
the RVars extension — load `RVars` alongside `DistributionPlots` to use it.
"""
function dimticks end

include("interface.jl")
include("density.jl")
include("intervals.jl")
include("slabs.jl")
include("dotlayout.jl")
include("geometry.jl")
include("recipes/slabinterval.jl")
include("recipes/dots.jl")
include("recipes/lineribbon.jl")

export asdist
export dimticks
export point_interval
export slab_curve
export dot_layout
export slabinterval, slabinterval!, SlabInterval
export halfeye, halfeye!, eye, eye!, ccdfinterval, ccdfinterval!, cdfinterval, cdfinterval!,
       gradientinterval, gradientinterval!, histinterval, histinterval!,
       slab, slab!, interval, interval!, pointinterval, pointinterval!, spike, spike!,
       HalfEye, Eye, CcdfInterval, CdfInterval, GradientInterval, HistInterval,
       Slab, Interval, PointInterval, Spike
export dots, dots!, dotsinterval, dotsinterval!, Dots, DotsInterval
export lineribbon, lineribbon!, LineRibbon

end # module
