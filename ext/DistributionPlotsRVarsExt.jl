module DistributionPlotsRVarsExt

using DistributionPlots
using DistributionPlots: _to_dist_args, AbstractDist, asdist
using RVars
using RVars: RVar, draws, variables
using Makie
using Distributions

# scalar RV (N=0): one distribution from all draws.
# vector RV (N=1, length k): k distributions at integer positions 1:k.
function DistributionPlots._to_dist_args(x::RVar)
    raw = draws(x)                       # (ndraws, dims...) — axis 1 is draws
    if ndims(raw) == 1
        return ([1.0], AbstractDist[asdist(vec(raw))])
    else
        k = size(raw, 2)
        return (collect(1.0:k), AbstractDist[asdist(vec(raw[:, j])) for j in 1:k])
    end
end

# Register convert_arguments for every public recipe type.
const _RECIPES = (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)

for T in _RECIPES
    @eval Makie.convert_arguments(::Type{<:$T}, x::RVar) = (_to_dist_args(x),)
end

end # module
