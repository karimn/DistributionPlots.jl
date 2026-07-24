module DistributionPlotsMCMCChainsExt

using DistributionPlots
using DistributionPlots: _to_dist_args
using MCMCChains: Chains
using RVars: RVar
using Makie

# Route Chains through RVars' own Chains conversion (single tested path),
# then reuse the RVars _to_dist_args method.
DistributionPlots._to_dist_args(chn::Chains) = _to_dist_args(RVar(chn))

const _RECIPES = (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)

for T in _RECIPES
    @eval Makie.convert_arguments(::Type{<:$T}, chn::Chains) = (_to_dist_args(chn),)
end

end # module
